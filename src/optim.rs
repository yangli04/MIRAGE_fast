//! Numerical optimization routines used by MIRAGE MLEs.
//!
//! - `brent_min`: Brent's algorithm for 1-D bounded minimization. Mirrors R's
//!   `optim(method="Brent", lower, upper)`.
//! - `lbfgsb_2d`: bounded 2-D minimization for the joint (beta, kappa) MLE.
//!   We implement a rectangular Nelder-Mead with bound clamping. R's
//!   L-BFGS-B converges to the same MLE because the joint likelihood for
//!   the heter model is well-behaved (smooth, finite, near-quadratic at the
//!   optimum). Outputs are rounded to 5 decimals to match the R reference.

/// Brent's method for univariate minimization on [lo, hi].
///
/// Translated from Numerical Recipes / R's optim Brent driver.
/// `tol` defaults follow R's optim: relative `1.49e-8`.
pub fn brent_min<F: FnMut(f64) -> f64>(mut f: F, lo: f64, hi: f64, tol: f64) -> f64 {
    let cgold = 0.3819660112501051_f64;
    let zeps = 1.0e-10_f64;
    let mut a = lo;
    let mut b = hi;
    let mut x = a + cgold * (b - a);
    let mut w = x;
    let mut v = x;
    let mut fx = f(x);
    let mut fw = fx;
    let mut fv = fx;
    let mut d = 0.0_f64;
    let mut e = 0.0_f64;

    for _ in 0..200 {
        let xm = 0.5 * (a + b);
        let tol1 = tol * x.abs() + zeps;
        let tol2 = 2.0 * tol1;
        if (x - xm).abs() <= (tol2 - 0.5 * (b - a)) {
            return x;
        }
        let mut use_golden = true;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let q0 = (x - v) * (fx - fw);
            let mut p = (x - v) * q0 - (x - w) * r;
            let mut q = 2.0 * (q0 - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let etemp = e;
            e = d;
            if p.abs() < (0.5 * q * etemp).abs() && p > q * (a - x) && p < q * (b - x) {
                d = p / q;
                let u = x + d;
                if (u - a) < tol2 || (b - u) < tol2 {
                    d = if xm - x >= 0.0 { tol1 } else { -tol1 };
                }
                use_golden = false;
            }
        }
        if use_golden {
            e = if x >= xm { a - x } else { b - x };
            d = cgold * e;
        }
        let u = if d.abs() >= tol1 {
            x + d
        } else if d >= 0.0 {
            x + tol1
        } else {
            x - tol1
        };
        let fu = f(u);
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    x
}

/// Nelder-Mead simplex with hard box-constraint clamping. Kept as a generic
/// 2D optimizer; the joint (beta, kappa) MLE uses coordinate descent with
/// Brent (matches L-BFGS-B more reliably than NM in this problem).
#[allow(dead_code)]
pub fn nelder_mead_2d<F: FnMut(f64, f64) -> f64>(
    mut f: F,
    init: (f64, f64),
    lo: (f64, f64),
    hi: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> (f64, f64) {
    let clamp = |p: (f64, f64)| -> (f64, f64) {
        (p.0.clamp(lo.0, hi.0), p.1.clamp(lo.1, hi.1))
    };
    let init = clamp(init);
    let step = (
        (0.05_f64).min((hi.0 - lo.0) * 0.1).max(1e-4),
        (0.05_f64).min((hi.1 - lo.1) * 0.1).max(1e-4),
    );
    let mut simplex: [(f64, f64); 3] = [
        init,
        clamp((init.0 + step.0, init.1)),
        clamp((init.0, init.1 + step.1)),
    ];
    let mut fvals: [f64; 3] = [
        f(simplex[0].0, simplex[0].1),
        f(simplex[1].0, simplex[1].1),
        f(simplex[2].0, simplex[2].1),
    ];

    let alpha = 1.0;
    let gamma = 2.0;
    let rho = 0.5;
    let sigma = 0.5;

    for _ in 0..max_iter {
        // Sort indices by fvals ascending.
        let mut idx = [0, 1, 2];
        idx.sort_by(|&a, &b| fvals[a].partial_cmp(&fvals[b]).unwrap_or(std::cmp::Ordering::Equal));
        let (best, mid, worst) = (idx[0], idx[1], idx[2]);
        let max_dist = (simplex[best].0 - simplex[worst].0)
            .abs()
            .max((simplex[best].1 - simplex[worst].1).abs());
        if max_dist < tol && (fvals[worst] - fvals[best]).abs() < tol {
            break;
        }
        // Centroid of best+mid.
        let c = (
            0.5 * (simplex[best].0 + simplex[mid].0),
            0.5 * (simplex[best].1 + simplex[mid].1),
        );
        // Reflection.
        let xr = clamp((c.0 + alpha * (c.0 - simplex[worst].0), c.1 + alpha * (c.1 - simplex[worst].1)));
        let fr = f(xr.0, xr.1);
        if fr < fvals[best] {
            // Expansion.
            let xe = clamp((c.0 + gamma * (xr.0 - c.0), c.1 + gamma * (xr.1 - c.1)));
            let fe = f(xe.0, xe.1);
            if fe < fr {
                simplex[worst] = xe;
                fvals[worst] = fe;
            } else {
                simplex[worst] = xr;
                fvals[worst] = fr;
            }
            continue;
        }
        if fr < fvals[mid] {
            simplex[worst] = xr;
            fvals[worst] = fr;
            continue;
        }
        // Contraction.
        let xc = clamp((c.0 + rho * (simplex[worst].0 - c.0), c.1 + rho * (simplex[worst].1 - c.1)));
        let fc = f(xc.0, xc.1);
        if fc < fvals[worst] {
            simplex[worst] = xc;
            fvals[worst] = fc;
            continue;
        }
        // Shrink towards best.
        for j in [mid, worst] {
            simplex[j] = clamp((
                simplex[best].0 + sigma * (simplex[j].0 - simplex[best].0),
                simplex[best].1 + sigma * (simplex[j].1 - simplex[best].1),
            ));
            fvals[j] = f(simplex[j].0, simplex[j].1);
        }
    }

    // Return the best vertex.
    let mut idx = [0, 1, 2];
    idx.sort_by(|&a, &b| fvals[a].partial_cmp(&fvals[b]).unwrap_or(std::cmp::Ordering::Equal));
    simplex[idx[0]]
}

/// Compute the median of a slice. NaNs are filtered out. Linear-time
/// nth_element-style would suffice but we just sort for clarity.
pub fn median(values: &[f64]) -> f64 {
    let mut v: Vec<f64> = values.iter().cloned().filter(|x| !x.is_nan()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Compute a quantile using R's default type-7 (continuous) interpolation,
/// matching `quantile(x, probs)` in R.
pub fn quantile_type7(values: &[f64], q: f64) -> f64 {
    let mut v: Vec<f64> = values.iter().cloned().filter(|x| !x.is_nan()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 1 {
        return v[0];
    }
    let h = (n as f64 - 1.0) * q;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = h - lo as f64;
    v[lo] * (1.0 - frac) + v[hi] * frac
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brent_min_quadratic() {
        // Minimum of (x - 0.3)^2 on [0,1] is 0.3.
        let x = brent_min(|x| (x - 0.3).powi(2), 0.0, 1.0, 1e-9);
        assert!((x - 0.3).abs() < 1e-6);
    }

    #[test]
    fn nm_quadratic() {
        // f(x,y) = (x-0.2)^2 + (y-0.7)^2
        let p = nelder_mead_2d(
            |x, y| (x - 0.2).powi(2) + (y - 0.7).powi(2),
            (0.5, 0.5),
            (1e-6, 1e-6),
            (1.0 - 1e-6, 1.0 - 1e-6),
            1e-9,
            500,
        );
        assert!((p.0 - 0.2).abs() < 1e-3, "x={}", p.0);
        assert!((p.1 - 0.7).abs() < 1e-3, "y={}", p.1);
    }

    #[test]
    fn quantile_matches_r() {
        // R: quantile(c(1,2,3,4,5), 0.95) = 4.8
        let q = quantile_type7(&[1.0, 2.0, 3.0, 4.0, 5.0], 0.95);
        assert!((q - 4.8).abs() < 1e-9);
    }
}
