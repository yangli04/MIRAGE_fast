//! Statistical primitives used by MIRAGE: binomial / Fisher / LRT tests
//! and Benjamini-Hochberg FDR adjustment.
//!
//! All tests mirror R's `binom.test(alternative="greater")`,
//! `fisher.test(alternative="greater")`, and the methylSig-style
//! log-likelihood ratio test used in the original R package.

use statrs::distribution::{Binomial, ChiSquared, ContinuousCDF, DiscreteCDF};
#[allow(unused_imports)]
use statrs::distribution::Discrete;

/// Survival function P(X >= k) for X ~ Binomial(n, p).
/// Equivalent to R's `pbinom(k - 1, n, p, lower.tail = FALSE)`.
pub fn binom_sf_ge(k: u64, n: u64, p: f64) -> f64 {
    if n == 0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    let p = p.clamp(0.0, 1.0);
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let dist = Binomial::new(p, n).expect("valid binomial");
    1.0 - dist.cdf(k - 1)
}

/// CDF P(X <= k) for X ~ Binomial(n, p). Matches R's `pbinom(k, n, p)`
/// where k may be any real (truncated to floor by statrs).
pub fn binom_cdf(k: f64, n: u64, p: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let p = p.clamp(0.0, 1.0);
    if p <= 0.0 {
        return 1.0;
    }
    if p >= 1.0 {
        return if k >= n as f64 { 1.0 } else { 0.0 };
    }
    if k < 0.0 {
        return 0.0;
    }
    let k = k.floor();
    if k >= n as f64 {
        return 1.0;
    }
    let k_u: u64 = k as u64;
    let dist = Binomial::new(p, n).expect("valid binomial");
    dist.cdf(k_u)
}

/// PMF P(X = k) for X ~ Binomial(n, p). Matches R's `dbinom(k, n, p)`.
pub fn binom_pmf(k: u64, n: u64, p: f64) -> f64 {
    if k > n {
        return 0.0;
    }
    if n == 0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    let p = p.clamp(0.0, 1.0);
    if p == 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if p == 1.0 {
        return if k == n { 1.0 } else { 0.0 };
    }
    // Use log-space for numerical stability.
    let ln_pmf =
        log_binom_coeff(n, k) + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln();
    ln_pmf.exp()
}

/// log of binomial coefficient C(n, k) using log-gamma.
fn log_binom_coeff(n: u64, k: u64) -> f64 {
    log_gamma((n + 1) as f64) - log_gamma((k + 1) as f64) - log_gamma((n - k + 1) as f64)
}

fn log_gamma(x: f64) -> f64 {
    statrs::function::gamma::ln_gamma(x)
}

/// One-sided binomial test, alternative = "greater".
/// p-value is P(X >= k | n, p).
pub fn binom_test_greater(k: u64, n: u64, p: f64) -> f64 {
    binom_sf_ge(k, n, p).clamp(0.0, 1.0)
}

/// Fisher's exact test on a 2x2 table, one-sided "greater".
///
/// R's `fisher.test(matrix(c(a,c,b,d), 2, 2), alternative="greater")` tests
/// whether the odds ratio for the (1,1) cell exceeds 1. The hypergeometric
/// p-value is P(X >= a) given the marginals.
///
/// Convention: the table passed in is
///   [a b]   = [treatment_X  treatment_fixed]
///   [c d]     [control_X    control_fixed]
/// (matching the R call `cbind(t_X, t_fix, c_X, c_fix)` then `matrix(x, 2, 2)`
/// which fills column-major: first column = (a, b), second = (c, d) with a=tX,
/// b=tFix, c=cX, d=cFix.) Here `a` is the test cell.
pub fn fisher_test_greater(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }
    let row1 = a + b; // sum of column with test cell row (k draws)
    let col1 = a + c; // total successes in population (K)

    // Hypergeometric: drawing `row1` items from population `n`,
    // of which `col1` are successes; X = a follows Hypergeometric(n, col1, row1).
    // R's convention for fisher.test "greater" computes P(X >= a).
    let max_k = row1.min(col1);
    let min_k = if row1 + col1 > n { row1 + col1 - n } else { 0 };

    if a > max_k {
        return 0.0;
    }
    if a <= min_k {
        return 1.0;
    }

    // Numerically: sum hypergeom PMF for x in a..=max_k.
    let mut log_pmfs: Vec<f64> = Vec::with_capacity((max_k - min_k + 1) as usize);
    for x in min_k..=max_k {
        // PMF(x) = C(col1, x) * C(n-col1, row1-x) / C(n, row1)
        let lp = log_binom_coeff(col1, x)
            + log_binom_coeff(n - col1, row1 - x)
            - log_binom_coeff(n, row1);
        log_pmfs.push(lp);
    }
    let max_lp = log_pmfs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = log_pmfs.iter().map(|&lp| (lp - max_lp).exp()).sum();
    let geq_idx = (a - min_k) as usize;
    let upper: f64 = log_pmfs[geq_idx..]
        .iter()
        .map(|&lp| (lp - max_lp).exp())
        .sum();
    (upper / total).clamp(0.0, 1.0)
}

/// Log-likelihood ratio test for a 2x2 binomial table, ported from
/// methylSig's `diff_binomial`. Returns (ratio, p-value).
///
/// The R implementation:
///   ratio = 2 * [c_X log(c_X/c_n) + c_F log(c_F/c_n)
///                + t_X log(t_X/t_n) + t_F log(t_F/t_n)
///                - (c_X+t_X) log((c_X+t_X)/(c_n+t_n))
///                - (c_F+t_F) log((c_F+t_F)/(c_n+t_n))]
/// p-value = P(chi2_1 >= ratio).
pub fn log_lik_ratio_test(
    treatment_total: u64,
    treatment_fixed: u64,
    control_total: u64,
    control_fixed: u64,
) -> (f64, f64) {
    let eps = 1e-100_f64;
    let cx = (control_total - control_fixed) as f64;
    let cf = control_fixed as f64;
    let tx = (treatment_total - treatment_fixed) as f64;
    let tf = treatment_fixed as f64;
    let cn = control_total as f64;
    let tn = treatment_total as f64;

    let ratio = 2.0
        * (cx * ((cx / cn) + eps).ln()
            + cf * ((cf / cn) + eps).ln()
            + tx * ((tx / tn) + eps).ln()
            + tf * ((tf / tn) + eps).ln()
            - (cx + tx) * (((cx + tx) / (cn + tn)) + eps).ln()
            - (cf + tf) * (((cf + tf) / (cn + tn)) + eps).ln());
    let p = if ratio <= 0.0 || !ratio.is_finite() {
        1.0
    } else {
        let chi = ChiSquared::new(1.0).expect("df=1");
        1.0 - chi.cdf(ratio)
    };
    (ratio, p.clamp(0.0, 1.0))
}

/// Benjamini-Hochberg FDR adjustment, matching R's
/// `p.adjust(p, method = "BH")`.
pub fn bh_adjust(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    if n == 0 {
        return Vec::new();
    }
    // Indices sorted by p ascending; NA preserved as-is.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        let pi = pvals[i];
        let pj = pvals[j];
        match (pi.is_nan(), pj.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => pi.partial_cmp(&pj).unwrap(),
        }
    });
    let n_non_na = pvals.iter().filter(|p| !p.is_nan()).count();
    let mut adj = vec![f64::NAN; n];
    if n_non_na == 0 {
        return adj;
    }
    // q_i = p_(i) * n / i; then cumulative min from the largest rank backwards;
    // then clamp to <= 1.
    let mut q = vec![0.0f64; n_non_na];
    for (rank0, &idx) in order.iter().take(n_non_na).enumerate() {
        let i = (rank0 + 1) as f64;
        q[rank0] = pvals[idx] * (n_non_na as f64) / i;
    }
    // cumulative min from the end.
    for k in (0..n_non_na - 1).rev() {
        if q[k + 1] < q[k] {
            q[k] = q[k + 1];
        }
    }
    for (rank0, &idx) in order.iter().take(n_non_na).enumerate() {
        adj[idx] = q[rank0].min(1.0);
    }
    adj
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol || ((a - b).abs() / b.abs().max(1e-12) < tol)
    }

    #[test]
    fn binom_test_greater_matches_r() {
        // R: binom.test(7, 10, p=0.5, alternative="greater")$p.value = 0.171875
        let p = binom_test_greater(7, 10, 0.5);
        assert!(approx(p, 0.171875, 1e-6), "got {}", p);
    }

    #[test]
    fn fisher_greater_matches_r() {
        // R: fisher.test(matrix(c(8,2,1,5),2,2), alternative="greater")$p.value
        //  = 0.024475524475524 (P(X >= 8 | N=16, K=10, n=9), hypergeom).
        let p = fisher_test_greater(8, 2, 1, 5);
        assert!(approx(p, 0.0244755244755, 1e-9), "got {}", p);
    }

    #[test]
    fn lrt_chisq() {
        let (ratio, p) = log_lik_ratio_test(10, 5, 10, 9);
        assert!(ratio > 0.0);
        assert!(p > 0.0 && p < 1.0);
    }

    #[test]
    fn bh_matches_r() {
        // R: p.adjust(c(0.01, 0.04, 0.03, 0.5), "BH")
        // = 0.040 0.053... 0.053... 0.500
        let a = bh_adjust(&[0.01, 0.04, 0.03, 0.5]);
        let exp = [0.04, 0.053333333, 0.053333333, 0.5];
        for (x, y) in a.iter().zip(exp.iter()) {
            assert!(approx(*x, *y, 1e-6), "got {:?}", a);
        }
    }
}
