//! Maximum-likelihood estimators and Bayesian inference helpers.
//!
//! Direct ports of:
//!   * `MLE_for_lambda1.R`
//!   * `MLE_for_lambda2.R`
//!   * `MLE_for_beta_homo.R`
//!   * `MLE_joint_for_beta_kappa_heter.R`
//!
//! The R implementation parameterizes a per-site mutation rate as
//!   p_mut = lambda * (1 - delta/3) + (1 - lambda) * delta
//! where `delta` is the sequencing error rate and `lambda` is the
//! true on-target conversion rate.

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

use crate::optim::brent_min;
use crate::stats::{binom_cdf, binom_pmf};

/// MLE for the high-signal mutation rate lambda1, given fixed/total counts in
/// the candidate high-signal site set.
pub fn mle_for_lambda1(
    n_fixed: &[u64],
    n_total: &[u64],
    delta: f64,
    initial_lambda1: f64,
) -> f64 {
    if n_fixed.is_empty() {
        // R behaviour: optim Brent on a constant likelihood returns the
        // initial parameter unchanged.
        return initial_lambda1;
    }
    let aa = 1.0 - delta / 3.0;
    let neg_log_lik = |lambda1: f64| -> f64 {
        let bb = lambda1 * aa + (1.0 - lambda1) * delta;
        if bb <= 0.0 || bb >= 1.0 {
            return f64::INFINITY;
        }
        let log_bb = bb.ln();
        let log_1mbb = (1.0 - bb).ln();
        let mut nt = 0.0;
        for (&f, &t) in n_fixed.iter().zip(n_total.iter()) {
            let f = f as f64;
            let t = t as f64;
            nt += (t - f) * log_bb + f * log_1mbb;
        }
        -nt
    };
    brent_min(neg_log_lik, 0.0, 1.0, 1.49e-8)
}

/// MLE for the background mutation rate lambda2. Subsamples up to
/// `num_sites` rows (matches the R reference).
pub fn mle_for_lambda2(
    n_fixed: &[u64],
    n_total: &[u64],
    num_sites: usize,
    delta: f64,
    _initial_lambda2: f64,
    seed: Option<u64>,
) -> f64 {
    let total_len = n_fixed.len();
    let take = num_sites.min(total_len);
    let (use_fixed, use_total): (Vec<u64>, Vec<u64>) = if take == total_len {
        (n_fixed.to_vec(), n_total.to_vec())
    } else {
        let mut rng: ChaCha20Rng = match seed {
            Some(s) => ChaCha20Rng::seed_from_u64(s),
            None => ChaCha20Rng::from_entropy(),
        };
        let mut idx: Vec<usize> = (0..total_len).collect();
        idx.shuffle(&mut rng);
        idx.truncate(take);
        let f = idx.iter().map(|&i| n_fixed[i]).collect();
        let t = idx.iter().map(|&i| n_total[i]).collect();
        (f, t)
    };

    let aa = 1.0 - delta / 3.0;
    let neg_log_lik = |lambda2: f64| -> f64 {
        let bb = lambda2 * aa + (1.0 - lambda2) * delta;
        if bb <= 0.0 || bb >= 1.0 {
            return f64::INFINITY;
        }
        let log_bb = bb.ln();
        let log_1mbb = (1.0 - bb).ln();
        let mut nt = 0.0;
        for (&f, &t) in use_fixed.iter().zip(use_total.iter()) {
            let f = f as f64;
            let t = t as f64;
            nt += (t - f) * log_bb + f * log_1mbb;
        }
        -nt
    };
    brent_min(neg_log_lik, 0.0, 1.0, 1.49e-8)
}

/// MLE for site-level signal beta in homozygous sites, given lambda1 and
/// lambda2 estimates.
pub fn mle_for_beta_homo(
    n_fixed_treatment: u64,
    n_total_treatment: u64,
    lambda1: f64,
    lambda2: f64,
    delta: f64,
    _initial_beta: f64,
) -> f64 {
    let f = n_fixed_treatment as f64;
    let t = n_total_treatment as f64;
    let neg_log_lik = |beta: f64| -> f64 {
        let aa = lambda1 * beta + lambda2 * (1.0 - beta);
        let bb = aa * (1.0 - delta / 3.0) + (1.0 - aa) * delta;
        if bb <= 0.0 || bb >= 1.0 {
            return f64::INFINITY;
        }
        -((t - f) * bb.ln() + f * (1.0 - bb).ln())
    };
    brent_min(neg_log_lik, 1e-6, 1.0, 1.49e-8)
}

/// Posterior probability and beta MLE for a homozygous site given a motif
/// prior. Returns (beta_est, posterior).
///
/// The "posterior" here is the posterior probability that the site is
/// **un-methylated**, matching the R convention:
///   posterior = (1 - prior) * Pr(D | unmethyl)
///             / [(1 - prior) * Pr(D | unmethyl) + prior * Pr(D | methyl)]
pub fn bayesian_inference_beta_homo(
    n_fixed_treatment: u64,
    n_total_treatment: u64,
    lambda1: f64,
    lambda2: f64,
    motif_prior: f64,
    delta: f64,
    initial_beta: f64,
) -> (f64, f64) {
    let aa = lambda2;
    let bb = aa * (1.0 - delta / 3.0) + (1.0 - aa) * delta;
    let unmethyl_likelihood = binom_pmf(
        n_total_treatment.saturating_sub(n_fixed_treatment),
        n_total_treatment,
        bb,
    );

    // homo_Likelihood from R: P(N_X | methyl, integrated over beta)
    let n_total = n_total_treatment as f64;
    let n_fix = n_fixed_treatment as f64;
    let p1_arg = (n_total + 1.0 - n_fix).max(0.0);
    let p1 = binom_cdf(p1_arg, n_total_treatment + 1, (1.0 - delta * 4.0 / 3.0) * lambda2)
        - binom_cdf(p1_arg, n_total_treatment + 1, (1.0 - delta * 4.0 / 3.0) * lambda1);
    let normalizing = (1.0 - delta * 4.0 / 3.0) * (lambda1 - lambda2) * (n_total + 1.0);
    let methyl_likelihood = if normalizing == 0.0 { 0.0 } else { p1 / normalizing };

    let a_term = (1.0 - motif_prior) * unmethyl_likelihood;
    let b_term = motif_prior * methyl_likelihood;
    let denom = a_term + b_term;
    let posterior = if denom == 0.0 || !denom.is_finite() {
        f64::NAN
    } else {
        a_term / denom
    };

    let beta = mle_for_beta_homo(
        n_fixed_treatment,
        n_total_treatment,
        lambda1,
        lambda2,
        delta,
        initial_beta,
    );
    (beta, posterior)
}

/// Joint MLE for (beta, kappa) on heterozygous sites.
/// Returns the rounded-to-5-decimal estimates, matching the R reference.
pub fn mle_joint_beta_kappa_heter(
    n_fixed_treatment: u64,
    n_total_treatment: u64,
    n_fixed_control: u64,
    n_total_control: u64,
    lambda1: f64,
    lambda2: f64,
    delta: f64,
    initial_beta: f64,
    initial_kappa: f64,
) -> (f64, f64) {
    let tf = n_fixed_treatment as f64;
    let tn = n_total_treatment as f64;
    let cf = n_fixed_control as f64;
    let cn = n_total_control as f64;
    let neg_log_lik = |beta: f64, kappa: f64| -> f64 {
        let aa = lambda1 * beta + lambda2 * (1.0 - beta);
        let bb = aa * (1.0 - delta / 3.0) + (1.0 - aa) * delta;
        let cc = (1.0 - kappa) * (1.0 - delta / 3.0);
        let p_t = bb * kappa + cc;
        let p_c = delta * kappa + cc;
        if p_t <= 0.0 || p_t >= 1.0 || p_c <= 0.0 || p_c >= 1.0 {
            return f64::INFINITY;
        }
        let nt = (tn - tf) * p_t.ln() + tf * (1.0 - p_t).ln()
            + (cn - cf) * p_c.ln() + cf * (1.0 - p_c).ln();
        -nt
    };

    let lo = 1e-6_f64;
    let hi = 1.0 - 1e-6_f64;
    let mut beta = initial_beta.clamp(lo, hi);
    let mut kappa = initial_kappa.clamp(lo, hi);
    let mut prev = f64::INFINITY;
    for _ in 0..200 {
        let neg_for_beta = |b: f64| neg_log_lik(b, kappa);
        beta = brent_min(neg_for_beta, lo, hi, 1.49e-9);
        let neg_for_kappa = |k: f64| neg_log_lik(beta, k);
        kappa = brent_min(neg_for_kappa, lo, hi, 1.49e-9);
        let cur = neg_log_lik(beta, kappa);
        if (prev - cur).abs() < 1e-12 {
            break;
        }
        prev = cur;
    }
    (round5(beta), round5(kappa))
}

/// MLE for kappa on heter sites using the control sample only.
pub fn mle_for_kappa_heter(
    n_fixed_control: u64,
    n_total_control: u64,
    delta: f64,
    _initial_kappa: f64,
) -> f64 {
    let f = n_fixed_control as f64;
    let t = n_total_control as f64;
    let neg_log_lik = |kappa: f64| -> f64 {
        let bb = (1.0 - kappa) * (1.0 - delta / 3.0) + kappa * delta;
        if bb <= 0.0 || bb >= 1.0 {
            return f64::INFINITY;
        }
        -((t - f) * bb.ln() + f * (1.0 - bb).ln())
    };
    brent_min(neg_log_lik, 1e-6, 1.0 - 1e-6, 1.49e-8)
}

/// MLE for beta on heter sites with kappa held fixed.
pub fn mle_for_beta_heter(
    n_fixed_treatment: u64,
    n_total_treatment: u64,
    kappa: f64,
    lambda1: f64,
    lambda2: f64,
    delta: f64,
    _initial_beta: f64,
) -> f64 {
    let f = n_fixed_treatment as f64;
    let t = n_total_treatment as f64;
    let neg_log_lik = |beta: f64| -> f64 {
        let aa = lambda1 * beta + lambda2 * (1.0 - beta);
        let bb = aa * (1.0 - delta / 3.0) + (1.0 - aa) * delta;
        let cc = (1.0 - kappa) * (1.0 - delta / 3.0);
        let dd = bb * kappa + cc;
        if dd <= 0.0 || dd >= 1.0 {
            return f64::INFINITY;
        }
        -((t - f) * dd.ln() + f * (1.0 - dd).ln())
    };
    brent_min(neg_log_lik, 1e-6, 1.0, 1.49e-8)
}

/// Bayesian inference for (kappa, beta) on heter sites with motif prior.
/// Returns (kappa_est, beta_est, posterior).
pub fn bayesian_inference_beta_kappa_heter(
    n_fixed_treatment: u64,
    n_total_treatment: u64,
    n_fixed_control: u64,
    n_total_control: u64,
    lambda1: f64,
    lambda2: f64,
    motif_prior: f64,
    delta: f64,
) -> (f64, f64, f64) {
    let kappa_est = mle_for_kappa_heter(n_fixed_control, n_total_control, delta, 0.5);

    let bb = lambda2 * (1.0 - delta / 3.0) + (1.0 - lambda2) * delta;
    let cc = (1.0 - kappa_est) * (1.0 - delta / 3.0);
    let dd = bb * kappa_est + cc;
    let ee = lambda1 * (1.0 - delta / 3.0) + (1.0 - lambda1) * delta;

    let unmethyl_likelihood = binom_pmf(
        n_total_treatment.saturating_sub(n_fixed_treatment),
        n_total_treatment,
        dd,
    );

    let n_total = n_total_treatment as f64;
    let n_fix = n_fixed_treatment as f64;
    let p2_arg = (n_total + 1.0 - n_fix).max(0.0);
    let prob2 = binom_cdf(p2_arg, n_total_treatment + 1, dd)
        - binom_cdf(p2_arg, n_total_treatment + 1, cc + ee * kappa_est);
    let normalizing = (1.0 - delta * 4.0 / 3.0) * (lambda1 - lambda2) * (n_total + 1.0);
    let methyl_likelihood = if normalizing == 0.0 { 0.0 } else { prob2 / normalizing };

    let a_term = (1.0 - motif_prior) * unmethyl_likelihood;
    let b_term = motif_prior * methyl_likelihood;
    let denom = a_term + b_term;
    let posterior = if denom == 0.0 || !denom.is_finite() {
        f64::NAN
    } else {
        a_term / denom
    };

    let beta_est = mle_for_beta_heter(
        n_fixed_treatment,
        n_total_treatment,
        kappa_est,
        lambda1,
        lambda2,
        delta,
        0.5,
    );
    (kappa_est, beta_est, posterior)
}

fn round5(x: f64) -> f64 {
    (x * 1e5).round() / 1e5
}
