//! High-level inference pipelines: `estimate_inference_with_empirical` and
//! `estimate_inference_with_prior`.

use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Result};
use rayon::prelude::*;

use crate::mle::{
    bayesian_inference_beta_homo, bayesian_inference_beta_kappa_heter, mle_for_beta_homo,
    mle_for_lambda1, mle_for_lambda2, mle_joint_beta_kappa_heter,
};
use crate::motif::generate_motif_sequences;
use crate::optim::{median, quantile_type7};
use crate::prior::MotifPriorRow;
use crate::stats::{
    bh_adjust, binom_test_greater, fisher_test_greater, log_lik_ratio_test,
};

/// One row of the input count table.
#[derive(Debug, Clone)]
pub struct SiteRecord {
    pub pos: String,
    pub motif: String,
    pub motif_type: String,
    pub treated_fixed_count: u64,
    pub treated_depth: u64,
    pub control_fixed_count: u64,
    pub control_depth: u64,
}

#[derive(Debug, Clone, Copy)]
pub enum BgMethod {
    Binomial,
    Fisher,
    Lrt,
}

#[derive(Debug, Clone, Copy)]
pub enum BgTarget {
    Treatment,
    Control,
    Both,
}

#[derive(Debug, Clone)]
pub enum Lambda1Mode {
    Auto,
    Site(Vec<String>),
    Fixed(f64),
}

#[derive(Debug, Clone, Copy)]
pub enum NmerMode {
    Mer5,
    F4,
    L4,
}

#[derive(Debug, Clone)]
pub struct EmpiricalParams {
    pub delta: f64,
    pub depth_cutoff: u64,
    pub homo_cutoff: f64,
    pub lambda1: Lambda1Mode,
    pub bg_method: BgMethod,
    pub bg_target: BgTarget,
    pub highly_methyl_cutoff: f64,
    pub seed: Option<u64>,
    pub thread: usize,
}

impl Default for EmpiricalParams {
    fn default() -> Self {
        Self {
            delta: 0.001,
            depth_cutoff: 10,
            homo_cutoff: 0.99,
            lambda1: Lambda1Mode::Auto,
            bg_method: BgMethod::Fisher,
            bg_target: BgTarget::Treatment,
            highly_methyl_cutoff: 0.95,
            seed: None,
            thread: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PriorParams {
    pub motif_specific: bool,
    pub delta: f64,
    pub depth_cutoff: u64,
    pub homo_cutoff: f64,
    pub bg_method: BgMethod,
    pub highly_methyl_cutoff: f64,
    pub seed: Option<u64>,
    pub thread: usize,
    pub motif: String,
    pub nmer: NmerMode,
}

impl Default for PriorParams {
    fn default() -> Self {
        Self {
            motif_specific: true,
            delta: 0.001,
            depth_cutoff: 10,
            homo_cutoff: 0.99,
            bg_method: BgMethod::Fisher,
            highly_methyl_cutoff: 0.95,
            seed: None,
            thread: 1,
            motif: "DRACH".into(),
            nmer: NmerMode::Mer5,
        }
    }
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct HomoSite {
    pub pos: String,
    pub motif: String,
    pub motif_type: String,
    pub treatment_x_rate: f64,
    pub control_x_rate: f64,
    pub treatment_fixed_rate: f64,
    pub control_fixed_rate: f64,
    pub treatment_fixed_count: u64,
    pub control_fixed_count: u64,
    pub treatment_depth: u64,
    pub control_depth: u64,
    pub binom_p: f64,
    pub binom_fdr: f64,
    pub fisher_p: f64,
    pub fisher_fdr: f64,
    pub log_lik_ratio: f64,
    pub lrt_p: f64,
    pub lrt_fdr: f64,
    pub beta_est: f64,
    pub lambda1: Option<f64>,
    pub prior_methylated: Option<f64>,
    pub posterior: Option<f64>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct HeterSite {
    pub pos: String,
    pub motif: String,
    pub motif_type: String,
    pub treatment_x_rate: f64,
    pub control_x_rate: f64,
    pub treatment_fixed_rate: f64,
    pub control_fixed_rate: f64,
    pub treatment_fixed_count: u64,
    pub control_fixed_count: u64,
    pub treatment_depth: u64,
    pub control_depth: u64,
    pub beta_est: f64,
    pub kappa_est: f64,
    pub binom_p: f64,
    pub binom_fdr: f64,
    pub fisher_p: f64,
    pub fisher_fdr: f64,
    pub log_lik_ratio: f64,
    pub lrt_p: f64,
    pub lrt_fdr: f64,
    pub lambda1: Option<f64>,
    pub prior_methylated: Option<f64>,
    pub posterior: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct EmpiricalResult {
    pub homosites: Vec<HomoSite>,
    pub hetersites: Vec<HeterSite>,
    pub lambda1: f64,
    pub lambda2: f64,
}

#[derive(Debug, Clone)]
pub struct PriorResult {
    pub homosites: Vec<HomoSite>,
    pub hetersites: Vec<HeterSite>,
    pub lambda1_summary: Vec<MotifPriorRow>,
    pub lambda1_per_motif: HashMap<String, f64>,
    pub lambda2: f64,
}

struct AllSitesView {
    pos: Vec<String>,
    motif: Vec<String>,
    motif_type: Vec<String>,
    t_fixed: Vec<u64>,
    t_total: Vec<u64>,
    c_fixed: Vec<u64>,
    c_total: Vec<u64>,
    t_x_rate: Vec<f64>,
    c_x_rate: Vec<f64>,
    t_fixed_rate: Vec<f64>,
    c_fixed_rate: Vec<f64>,
}

fn build_allsites(records: &[SiteRecord], depth_cutoff: u64) -> AllSitesView {
    let mut a = AllSitesView {
        pos: vec![],
        motif: vec![],
        motif_type: vec![],
        t_fixed: vec![],
        t_total: vec![],
        c_fixed: vec![],
        c_total: vec![],
        t_x_rate: vec![],
        c_x_rate: vec![],
        t_fixed_rate: vec![],
        c_fixed_rate: vec![],
    };
    for r in records {
        if r.treated_depth < depth_cutoff || r.control_depth < depth_cutoff {
            continue;
        }
        let tfr = r.treated_fixed_count as f64 / r.treated_depth as f64;
        let cfr = r.control_fixed_count as f64 / r.control_depth as f64;
        a.pos.push(r.pos.clone());
        a.motif.push(r.motif.clone());
        a.motif_type.push(r.motif_type.clone());
        a.t_fixed.push(r.treated_fixed_count);
        a.t_total.push(r.treated_depth);
        a.c_fixed.push(r.control_fixed_count);
        a.c_total.push(r.control_depth);
        a.t_fixed_rate.push(tfr);
        a.c_fixed_rate.push(cfr);
        a.t_x_rate.push(1.0 - tfr);
        a.c_x_rate.push(1.0 - cfr);
    }
    a
}

#[derive(Debug, Clone)]
struct TestRow {
    binom_p: f64,
    binom_fdr: f64,
    fisher_p: f64,
    fisher_fdr: f64,
    lrt: f64,
    lrt_p: f64,
    lrt_fdr: f64,
}

fn test_for_mutation_given_a_rate(
    treatment_fixed: &[u64],
    treatment_total: &[u64],
    control_fixed: &[u64],
    control_total: &[u64],
    mutation_rate: f64,
    thread: usize,
) -> Vec<TestRow> {
    let n = treatment_fixed.len();
    let make_tests = |i: usize| -> (f64, f64, f64, f64) {
        let tn = treatment_total[i];
        let tf = treatment_fixed[i];
        let cn = control_total[i];
        let cf = control_fixed[i];
        let tx = tn - tf;
        let cx = cn - cf;
        let bp = binom_test_greater(tx, tn, mutation_rate);
        // fisher.test on matrix(c(tx, tf, cx, cf), 2, 2) "greater"
        let fp = fisher_test_greater(tx, tf, cx, cf);
        let (ratio, lp) = log_lik_ratio_test(tn, tf, cn, cf);
        (bp, fp, ratio, lp)
    };
    let stats: Vec<(f64, f64, f64, f64)> = if thread > 1 {
        (0..n).into_par_iter().map(make_tests).collect()
    } else {
        (0..n).map(make_tests).collect()
    };
    let binom_p: Vec<f64> = stats.iter().map(|s| s.0).collect();
    let fisher_p: Vec<f64> = stats.iter().map(|s| s.1).collect();
    let lrt_ratio: Vec<f64> = stats.iter().map(|s| s.2).collect();
    let lrt_p: Vec<f64> = stats.iter().map(|s| s.3).collect();
    let binom_fdr = bh_adjust(&binom_p);
    let fisher_fdr = bh_adjust(&fisher_p);
    let lrt_fdr = bh_adjust(&lrt_p);
    (0..n)
        .map(|i| TestRow {
            binom_p: binom_p[i],
            binom_fdr: binom_fdr[i],
            fisher_p: fisher_p[i],
            fisher_fdr: fisher_fdr[i],
            lrt: lrt_ratio[i],
            lrt_p: lrt_p[i],
            lrt_fdr: lrt_fdr[i],
        })
        .collect()
}

fn split_homo<'a>(allsites: &'a AllSitesView, homo_cutoff: f64) -> (Vec<usize>, Vec<usize>) {
    let mut homo: Vec<usize> = Vec::new();
    let mut heter: Vec<usize> = Vec::new();
    for i in 0..allsites.pos.len() {
        if allsites.c_fixed_rate[i] > homo_cutoff {
            homo.push(i);
        } else {
            heter.push(i);
        }
    }
    (homo, heter)
}

/// Reproduce R's `quantile(treatment_X_rate, 0.95)` then take the median
/// of values below or equal to that cutoff — used as the initial background
/// rate.
fn initial_background_rate(t_x_rates: &[f64]) -> f64 {
    let q = quantile_type7(t_x_rates, 0.95);
    let below: Vec<f64> = t_x_rates.iter().cloned().filter(|x| *x <= q).collect();
    median(&below)
}

fn pick_indices_by_method(
    summary: &[TestRow],
    method: BgMethod,
) -> (Vec<usize>, Vec<usize>) {
    let mut g0 = Vec::new();
    let mut g1 = Vec::new();
    for (i, r) in summary.iter().enumerate() {
        let (p, fdr) = match method {
            BgMethod::Binomial => (r.binom_p, r.binom_fdr),
            BgMethod::Fisher => (r.fisher_p, r.fisher_fdr),
            BgMethod::Lrt => (r.lrt_p, r.lrt_fdr),
        };
        if !p.is_nan() && p > 0.05 {
            g0.push(i);
        }
        if !fdr.is_nan() && fdr < 0.05 {
            g1.push(i);
        }
    }
    (g0, g1)
}

/// Public: empirical inference (single rate per dataset).
pub fn estimate_inference_with_empirical(
    records: &[SiteRecord],
    params: &EmpiricalParams,
) -> Result<EmpiricalResult> {
    let allsites = build_allsites(records, params.depth_cutoff);
    let (homo_idx, heter_idx) = split_homo(&allsites, params.homo_cutoff);

    eprintln!(
        "Considering {} ({:.2}%) sites as homozygous...",
        homo_idx.len(),
        if !allsites.pos.is_empty() {
            100.0 * homo_idx.len() as f64 / allsites.pos.len() as f64
        } else {
            0.0
        }
    );

    let homo_t_x_rate: Vec<f64> = homo_idx.iter().map(|&i| allsites.t_x_rate[i]).collect();
    let homo_t_fixed: Vec<u64> = homo_idx.iter().map(|&i| allsites.t_fixed[i]).collect();
    let homo_t_total: Vec<u64> = homo_idx.iter().map(|&i| allsites.t_total[i]).collect();
    let homo_c_fixed: Vec<u64> = homo_idx.iter().map(|&i| allsites.c_fixed[i]).collect();
    let homo_c_total: Vec<u64> = homo_idx.iter().map(|&i| allsites.c_total[i]).collect();

    let initial_bg = initial_background_rate(&homo_t_x_rate);
    let initial_summary = test_for_mutation_given_a_rate(
        &homo_t_fixed,
        &homo_t_total,
        &homo_c_fixed,
        &homo_c_total,
        initial_bg,
        params.thread,
    );

    let (igroup0, igroup1) = pick_indices_by_method(&initial_summary, params.bg_method);
    eprintln!(
        "Considering {} ({:.2}%) sites as background (homozygous sites only)...",
        igroup0.len(),
        if !homo_idx.is_empty() {
            100.0 * igroup0.len() as f64 / homo_idx.len() as f64
        } else {
            0.0
        }
    );

    // Estimate lambda2.
    let (n_fixed_bg, n_total_bg): (Vec<u64>, Vec<u64>) = match params.bg_target {
        BgTarget::Treatment => {
            let f: Vec<u64> = igroup0.iter().map(|&i| homo_t_fixed[i]).collect();
            let t: Vec<u64> = igroup0.iter().map(|&i| homo_t_total[i]).collect();
            (f, t)
        }
        BgTarget::Control => {
            let f: Vec<u64> = igroup0.iter().map(|&i| homo_c_fixed[i]).collect();
            let t: Vec<u64> = igroup0.iter().map(|&i| homo_c_total[i]).collect();
            (f, t)
        }
        BgTarget::Both => {
            let mut f: Vec<u64> = igroup0.iter().map(|&i| homo_t_fixed[i]).collect();
            let mut t: Vec<u64> = igroup0.iter().map(|&i| homo_t_total[i]).collect();
            f.extend(igroup0.iter().map(|&i| homo_c_fixed[i]));
            t.extend(igroup0.iter().map(|&i| homo_c_total[i]));
            (f, t)
        }
    };

    let lambda2 = mle_for_lambda2(
        &n_fixed_bg,
        &n_total_bg,
        100_000,
        params.delta,
        0.1,
        params.seed,
    );

    let final_bg = lambda2 * (1.0 - params.delta / 3.0) + (1.0 - lambda2) * params.delta;

    let final_summary = test_for_mutation_given_a_rate(
        &homo_t_fixed,
        &homo_t_total,
        &homo_c_fixed,
        &homo_c_total,
        final_bg,
        params.thread,
    );

    // Estimate lambda1.
    let lambda1 = match &params.lambda1 {
        Lambda1Mode::Fixed(v) => *v,
        Lambda1Mode::Auto => {
            let high_idx: Vec<usize> = igroup1.clone();
            eprintln!(
                "Considering {} ({:.2}%) sites as high-signal site candidates (homozygous sites only)...",
                high_idx.len(),
                if !homo_idx.is_empty() {
                    100.0 * high_idx.len() as f64 / homo_idx.len() as f64
                } else {
                    0.0
                }
            );
            let high_rates: Vec<f64> = high_idx.iter().map(|&i| homo_t_x_rate[i]).collect();
            let high_cutoff = quantile_type7(&high_rates, params.highly_methyl_cutoff);
            let kept: Vec<usize> = high_idx
                .iter()
                .copied()
                .filter(|&i| homo_t_x_rate[i] >= high_cutoff)
                .collect();
            let f: Vec<u64> = kept.iter().map(|&i| homo_t_fixed[i]).collect();
            let t: Vec<u64> = kept.iter().map(|&i| homo_t_total[i]).collect();
            mle_for_lambda1(&f, &t, params.delta, 0.2)
        }
        Lambda1Mode::Site(top_sites) => {
            let set: HashSet<&String> = top_sites.iter().collect();
            let kept: Vec<usize> = (0..allsites.pos.len())
                .filter(|&i| set.contains(&allsites.pos[i]))
                .collect();
            eprintln!(
                "Considering {} sites as high-signal site candidates...",
                kept.len()
            );
            let f: Vec<u64> = kept.iter().map(|&i| allsites.t_fixed[i]).collect();
            let t: Vec<u64> = kept.iter().map(|&i| allsites.t_total[i]).collect();
            mle_for_lambda1(&f, &t, params.delta, 0.2)
        }
    };

    // Per-site beta MLE — homo.
    let beta_homo: Vec<f64> = if params.thread > 1 {
        (0..homo_idx.len())
            .into_par_iter()
            .map(|j| {
                mle_for_beta_homo(
                    homo_t_fixed[j],
                    homo_t_total[j],
                    lambda1,
                    lambda2,
                    params.delta,
                    0.2,
                )
            })
            .collect()
    } else {
        (0..homo_idx.len())
            .map(|j| {
                mle_for_beta_homo(
                    homo_t_fixed[j],
                    homo_t_total[j],
                    lambda1,
                    lambda2,
                    params.delta,
                    0.2,
                )
            })
            .collect()
    };

    let mut homosites: Vec<HomoSite> = Vec::with_capacity(homo_idx.len());
    for (j, &i) in homo_idx.iter().enumerate() {
        let s = &final_summary[j];
        homosites.push(HomoSite {
            pos: allsites.pos[i].clone(),
            motif: allsites.motif[i].clone(),
            motif_type: allsites.motif_type[i].clone(),
            treatment_x_rate: allsites.t_x_rate[i],
            control_x_rate: allsites.c_x_rate[i],
            treatment_fixed_rate: allsites.t_fixed_rate[i],
            control_fixed_rate: allsites.c_fixed_rate[i],
            treatment_fixed_count: allsites.t_fixed[i],
            control_fixed_count: allsites.c_fixed[i],
            treatment_depth: allsites.t_total[i],
            control_depth: allsites.c_total[i],
            binom_p: s.binom_p,
            binom_fdr: s.binom_fdr,
            fisher_p: s.fisher_p,
            fisher_fdr: s.fisher_fdr,
            log_lik_ratio: s.lrt,
            lrt_p: s.lrt_p,
            lrt_fdr: s.lrt_fdr,
            beta_est: beta_homo[j],
            lambda1: None,
            prior_methylated: None,
            posterior: None,
        });
    }

    // Heter.
    let heter_t_fixed: Vec<u64> = heter_idx.iter().map(|&i| allsites.t_fixed[i]).collect();
    let heter_t_total: Vec<u64> = heter_idx.iter().map(|&i| allsites.t_total[i]).collect();
    let heter_c_fixed: Vec<u64> = heter_idx.iter().map(|&i| allsites.c_fixed[i]).collect();
    let heter_c_total: Vec<u64> = heter_idx.iter().map(|&i| allsites.c_total[i]).collect();

    let heter_estimates: Vec<(f64, f64)> = if params.thread > 1 {
        (0..heter_idx.len())
            .into_par_iter()
            .map(|j| {
                mle_joint_beta_kappa_heter(
                    heter_t_fixed[j],
                    heter_t_total[j],
                    heter_c_fixed[j],
                    heter_c_total[j],
                    lambda1,
                    lambda2,
                    params.delta,
                    0.2,
                    0.5,
                )
            })
            .collect()
    } else {
        (0..heter_idx.len())
            .map(|j| {
                mle_joint_beta_kappa_heter(
                    heter_t_fixed[j],
                    heter_t_total[j],
                    heter_c_fixed[j],
                    heter_c_total[j],
                    lambda1,
                    lambda2,
                    params.delta,
                    0.2,
                    0.5,
                )
            })
            .collect()
    };
    let heter_beta: Vec<f64> = heter_estimates.iter().map(|x| x.0).collect();
    let heter_kappa: Vec<f64> = heter_estimates.iter().map(|x| x.1).collect();

    let heter_bg_rate: Vec<f64> = (0..heter_idx.len())
        .map(|j| {
            let kappa = heter_kappa[j];
            (1.0 - kappa) * (1.0 - params.delta / 3.0)
                + kappa * (lambda2 * (1.0 - params.delta / 3.0) + (1.0 - lambda2) * params.delta)
        })
        .collect();

    // Per-site test against per-site rate.
    let n_het = heter_idx.len();
    let make_tests_het = |j: usize| -> (f64, f64, f64, f64) {
        let tn = heter_t_total[j];
        let tf = heter_t_fixed[j];
        let cn = heter_c_total[j];
        let cf = heter_c_fixed[j];
        let mr = heter_bg_rate[j];
        let tx = tn - tf;
        let cx = cn - cf;
        let bp = binom_test_greater(tx, tn, mr);
        let fp = fisher_test_greater(tx, tf, cx, cf);
        let (ratio, lp) = log_lik_ratio_test(tn, tf, cn, cf);
        (bp, fp, ratio, lp)
    };
    let stats_het: Vec<(f64, f64, f64, f64)> = if params.thread > 1 {
        (0..n_het).into_par_iter().map(make_tests_het).collect()
    } else {
        (0..n_het).map(make_tests_het).collect()
    };
    let binom_p_h: Vec<f64> = stats_het.iter().map(|s| s.0).collect();
    let fisher_p_h: Vec<f64> = stats_het.iter().map(|s| s.1).collect();
    let lrt_h: Vec<f64> = stats_het.iter().map(|s| s.2).collect();
    let lrt_p_h: Vec<f64> = stats_het.iter().map(|s| s.3).collect();
    let binom_fdr_h = bh_adjust(&binom_p_h);
    let fisher_fdr_h = bh_adjust(&fisher_p_h);
    let lrt_fdr_h = bh_adjust(&lrt_p_h);

    let mut hetersites: Vec<HeterSite> = Vec::with_capacity(heter_idx.len());
    for (j, &i) in heter_idx.iter().enumerate() {
        hetersites.push(HeterSite {
            pos: allsites.pos[i].clone(),
            motif: allsites.motif[i].clone(),
            motif_type: allsites.motif_type[i].clone(),
            treatment_x_rate: allsites.t_x_rate[i],
            control_x_rate: allsites.c_x_rate[i],
            treatment_fixed_rate: allsites.t_fixed_rate[i],
            control_fixed_rate: allsites.c_fixed_rate[i],
            treatment_fixed_count: allsites.t_fixed[i],
            control_fixed_count: allsites.c_fixed[i],
            treatment_depth: allsites.t_total[i],
            control_depth: allsites.c_total[i],
            beta_est: heter_beta[j],
            kappa_est: heter_kappa[j],
            binom_p: binom_p_h[j],
            binom_fdr: binom_fdr_h[j],
            fisher_p: fisher_p_h[j],
            fisher_fdr: fisher_fdr_h[j],
            log_lik_ratio: lrt_h[j],
            lrt_p: lrt_p_h[j],
            lrt_fdr: lrt_fdr_h[j],
            lambda1: None,
            prior_methylated: None,
            posterior: None,
        });
    }

    Ok(EmpiricalResult {
        homosites,
        hetersites,
        lambda1,
        lambda2,
    })
}

/// Public: prior-aware inference.
pub fn estimate_inference_with_prior(
    records: &[SiteRecord],
    params: &PriorParams,
    ref_freq_tab: &[MotifPriorRow],
) -> Result<PriorResult> {
    if ref_freq_tab.is_empty() {
        return Err(anyhow!("ref_freq_tab is empty"));
    }
    let allsites = build_allsites(records, params.depth_cutoff);
    let (homo_idx, heter_idx) = split_homo(&allsites, params.homo_cutoff);

    eprintln!(
        "Considering {} ({:.2}%) sites as homozygous...",
        homo_idx.len(),
        if !allsites.pos.is_empty() {
            100.0 * homo_idx.len() as f64 / allsites.pos.len() as f64
        } else {
            0.0
        }
    );

    let homo_t_x_rate: Vec<f64> = homo_idx.iter().map(|&i| allsites.t_x_rate[i]).collect();
    let homo_t_fixed: Vec<u64> = homo_idx.iter().map(|&i| allsites.t_fixed[i]).collect();
    let homo_t_total: Vec<u64> = homo_idx.iter().map(|&i| allsites.t_total[i]).collect();
    let homo_c_fixed: Vec<u64> = homo_idx.iter().map(|&i| allsites.c_fixed[i]).collect();
    let homo_c_total: Vec<u64> = homo_idx.iter().map(|&i| allsites.c_total[i]).collect();
    let homo_motif: Vec<String> = homo_idx.iter().map(|&i| allsites.motif[i].clone()).collect();

    let initial_bg = initial_background_rate(&homo_t_x_rate);
    let initial_summary = test_for_mutation_given_a_rate(
        &homo_t_fixed,
        &homo_t_total,
        &homo_c_fixed,
        &homo_c_total,
        initial_bg,
        params.thread,
    );

    let (igroup0, igroup1) = pick_indices_by_method(&initial_summary, params.bg_method);
    eprintln!(
        "Considering {} ({:.2}%) sites as background (homozygous sites only)...",
        igroup0.len(),
        if !homo_idx.is_empty() {
            100.0 * igroup0.len() as f64 / homo_idx.len() as f64
        } else {
            0.0
        }
    );

    // Estimate lambda2 using treatment counts only (matches R for prior version).
    let n_fixed_bg: Vec<u64> = igroup0.iter().map(|&i| homo_t_fixed[i]).collect();
    let n_total_bg: Vec<u64> = igroup0.iter().map(|&i| homo_t_total[i]).collect();
    let lambda2 = mle_for_lambda2(
        &n_fixed_bg,
        &n_total_bg,
        100_000,
        params.delta,
        0.1,
        params.seed,
    );

    // Overall lambda1 from quantile-cutoff highly-methylated set.
    let high_idx0: Vec<usize> = igroup1.clone();
    eprintln!(
        "Considering {} ({:.2}%) sites as high-signal site candidates (homozygous sites only)...",
        high_idx0.len(),
        if !homo_idx.is_empty() {
            100.0 * high_idx0.len() as f64 / homo_idx.len() as f64
        } else {
            0.0
        }
    );
    let high_rates: Vec<f64> = high_idx0.iter().map(|&i| homo_t_x_rate[i]).collect();
    let high_cutoff = quantile_type7(&high_rates, params.highly_methyl_cutoff);
    let high_kept: Vec<usize> = high_idx0
        .iter()
        .copied()
        .filter(|&i| homo_t_x_rate[i] >= high_cutoff)
        .collect();
    let f_overall: Vec<u64> = high_kept.iter().map(|&i| homo_t_fixed[i]).collect();
    let t_overall: Vec<u64> = high_kept.iter().map(|&i| homo_t_total[i]).collect();
    let lambda1_overall = mle_for_lambda1(&f_overall, &t_overall, params.delta, 0.2);

    // Per-motif lambda1.
    let target_seqs = generate_motif_sequences(&params.motif)?;
    let target_set: HashSet<&str> = target_seqs.iter().map(|s| s.as_str()).collect();
    let mut lambda1_per_motif: HashMap<String, f64> = HashMap::new();

    let group_for = |motif: &str| -> Option<String> {
        match params.nmer {
            NmerMode::Mer5 => {
                if target_set.contains(motif) {
                    Some(motif.to_string())
                } else {
                    None
                }
            }
            NmerMode::F4 => {
                if motif.len() < 4 {
                    return None;
                }
                let key: String = motif.chars().take(4).collect();
                let any_match = target_seqs.iter().any(|s| s.starts_with(&key));
                if any_match {
                    Some(key)
                } else {
                    None
                }
            }
            NmerMode::L4 => {
                if motif.len() < 5 {
                    return None;
                }
                let key: String = motif.chars().skip(1).take(4).collect();
                let any_match = target_seqs.iter().any(|s| s.ends_with(&key));
                if any_match {
                    Some(key)
                } else {
                    None
                }
            }
        }
    };

    if params.motif_specific {
        // Group high-signal sites by their grouping key.
        let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
        for &i in &high_idx0 {
            if let Some(k) = group_for(&homo_motif[i]) {
                groups.entry(k).or_default().push(i);
            }
        }
        let mut group_lambda: HashMap<String, f64> = HashMap::new();
        for (k, idxs) in &groups {
            let rates: Vec<f64> = idxs.iter().map(|&i| homo_t_x_rate[i]).collect();
            let cutoff = quantile_type7(&rates, params.highly_methyl_cutoff);
            let kept: Vec<usize> = idxs
                .iter()
                .copied()
                .filter(|&i| homo_t_x_rate[i] >= cutoff)
                .collect();
            let f: Vec<u64> = kept.iter().map(|&i| homo_t_fixed[i]).collect();
            let t: Vec<u64> = kept.iter().map(|&i| homo_t_total[i]).collect();
            let est = mle_for_lambda1(&f, &t, params.delta, 0.2);
            group_lambda.insert(k.clone(), est);
        }
        for seq in &target_seqs {
            let key = group_for(seq).unwrap_or_else(|| seq.clone());
            let est = group_lambda.get(&key).copied().unwrap_or(f64::NAN);
            lambda1_per_motif.insert(seq.clone(), est);
        }
    } else {
        for seq in &target_seqs {
            lambda1_per_motif.insert(seq.clone(), lambda1_overall);
        }
    }

    // Non-target motifs use the overall estimate.
    let observed_motifs: HashSet<String> = allsites.motif.iter().cloned().collect();
    for m in &observed_motifs {
        if !target_set.contains(m.as_str()) {
            lambda1_per_motif.insert(m.clone(), lambda1_overall);
        }
    }

    // Build lambda1.summary as a join with the prior reference table.
    let prior_index: HashMap<&str, &MotifPriorRow> =
        ref_freq_tab.iter().map(|r| (r.motif.as_str(), r)).collect();
    let mut summary: Vec<MotifPriorRow> = Vec::new();
    let mut summary_motifs: HashSet<String> = HashSet::new();
    for (m, &lam) in &lambda1_per_motif {
        let row = prior_index.get(m.as_str()).cloned();
        let (tx_freq, genome_freq, freq, motif_type, prior_methylated) = match row {
            Some(r) => (
                r.tx_freq,
                r.genome_freq,
                r.freq,
                r.motif_type.clone(),
                r.prior_methylated,
            ),
            None => (f64::NAN, f64::NAN, f64::NAN, String::new(), f64::NAN),
        };
        summary.push(MotifPriorRow {
            motif: m.clone(),
            tx_freq,
            genome_freq,
            freq,
            motif_type,
            prior_methylated,
        });
        summary.last_mut().unwrap().prior_methylated = prior_methylated;
        summary_motifs.insert(m.clone());
        // Repurpose freq for lambda1 storage isn't ideal; leave it as is but
        // remember lambda1 separately:
        let _ = lam;
    }
    // Keep a parallel lambda1 vector inline within the summary by storing in
    // an extra map output below.

    // Bayesian inference for homo sites.
    let homo_motif_ref: Vec<&str> = homo_motif.iter().map(|s| s.as_str()).collect();
    let homo_lambda1: Vec<Option<f64>> = homo_motif_ref
        .iter()
        .map(|m| lambda1_per_motif.get(*m).copied())
        .collect();
    let homo_prior: Vec<Option<f64>> = homo_motif_ref
        .iter()
        .map(|m| prior_index.get(*m).map(|r| r.prior_methylated))
        .collect();

    // Mirror R: keep rows whose prior_methylated is non-NA. Rows with NA
    // lambda1 stay but produce NaN beta/posterior.
    let mut homosites: Vec<HomoSite> = Vec::new();
    let homo_kept: Vec<usize> = (0..homo_idx.len())
        .filter(|&j| matches!(homo_prior[j], Some(p) if !p.is_nan()))
        .collect();

    let compute_one = |j: usize| -> (usize, f64, f64, f64, f64) {
        let lam1 = homo_lambda1[j].unwrap_or(f64::NAN);
        let prior = homo_prior[j].unwrap_or(f64::NAN);
        if lam1.is_nan() {
            return (j, f64::NAN, f64::NAN, lam1, prior);
        }
        let (beta, posterior) = bayesian_inference_beta_homo(
            homo_t_fixed[j],
            homo_t_total[j],
            lam1,
            lambda2,
            prior,
            params.delta,
            0.2,
        );
        (j, beta, posterior, lam1, prior)
    };
    let computed: Vec<(usize, f64, f64, f64, f64)> = if params.thread > 1 {
        homo_kept.par_iter().map(|&j| compute_one(j)).collect()
    } else {
        homo_kept.iter().map(|&j| compute_one(j)).collect()
    };

    // We need binom_p / fisher_p / lrt fields too. Reuse initial_summary for
    // those (R uses the unaltered homosites table without re-testing here).
    // To match R's behaviour the estimate_inference_with_prior function does
    // not run a final test on homo (only the first one). We fill stat fields
    // with NaN in that case to match: actually R doesn't store any stats here,
    // it overwrites homosites without those columns. We'll fill with NaN.
    for (j, beta, posterior, lam1, prior) in computed {
        let i = homo_idx[j];
        homosites.push(HomoSite {
            pos: allsites.pos[i].clone(),
            motif: allsites.motif[i].clone(),
            motif_type: allsites.motif_type[i].clone(),
            treatment_x_rate: allsites.t_x_rate[i],
            control_x_rate: allsites.c_x_rate[i],
            treatment_fixed_rate: allsites.t_fixed_rate[i],
            control_fixed_rate: allsites.c_fixed_rate[i],
            treatment_fixed_count: allsites.t_fixed[i],
            control_fixed_count: allsites.c_fixed[i],
            treatment_depth: allsites.t_total[i],
            control_depth: allsites.c_total[i],
            binom_p: f64::NAN,
            binom_fdr: f64::NAN,
            fisher_p: f64::NAN,
            fisher_fdr: f64::NAN,
            log_lik_ratio: f64::NAN,
            lrt_p: f64::NAN,
            lrt_fdr: f64::NAN,
            beta_est: beta,
            lambda1: Some(lam1),
            prior_methylated: Some(prior),
            posterior: Some(posterior),
        });
    }

    // Heter Bayesian inference.
    let heter_t_fixed: Vec<u64> = heter_idx.iter().map(|&i| allsites.t_fixed[i]).collect();
    let heter_t_total: Vec<u64> = heter_idx.iter().map(|&i| allsites.t_total[i]).collect();
    let heter_c_fixed: Vec<u64> = heter_idx.iter().map(|&i| allsites.c_fixed[i]).collect();
    let heter_c_total: Vec<u64> = heter_idx.iter().map(|&i| allsites.c_total[i]).collect();
    let heter_motif: Vec<String> = heter_idx.iter().map(|&i| allsites.motif[i].clone()).collect();

    let heter_lambda1: Vec<Option<f64>> = heter_motif
        .iter()
        .map(|m| lambda1_per_motif.get(m.as_str()).copied())
        .collect();
    let heter_prior: Vec<Option<f64>> = heter_motif
        .iter()
        .map(|m| prior_index.get(m.as_str()).map(|r| r.prior_methylated))
        .collect();

    let heter_estimates: Vec<(usize, f64, f64, f64, Option<f64>, Option<f64>)> = if params.thread > 1 {
        (0..heter_idx.len())
            .into_par_iter()
            .map(|j| {
                let lam1 = heter_lambda1[j].unwrap_or(f64::NAN);
                let prior = heter_prior[j].unwrap_or(f64::NAN);
                if lam1.is_nan() || prior.is_nan() {
                    (j, f64::NAN, f64::NAN, f64::NAN, heter_lambda1[j], heter_prior[j])
                } else {
                    let (kappa, beta, posterior) = bayesian_inference_beta_kappa_heter(
                        heter_t_fixed[j],
                        heter_t_total[j],
                        heter_c_fixed[j],
                        heter_c_total[j],
                        lam1,
                        lambda2,
                        prior,
                        params.delta,
                    );
                    (j, kappa, beta, posterior, Some(lam1), Some(prior))
                }
            })
            .collect()
    } else {
        (0..heter_idx.len())
            .map(|j| {
                let lam1 = heter_lambda1[j].unwrap_or(f64::NAN);
                let prior = heter_prior[j].unwrap_or(f64::NAN);
                if lam1.is_nan() || prior.is_nan() {
                    (j, f64::NAN, f64::NAN, f64::NAN, heter_lambda1[j], heter_prior[j])
                } else {
                    let (kappa, beta, posterior) = bayesian_inference_beta_kappa_heter(
                        heter_t_fixed[j],
                        heter_t_total[j],
                        heter_c_fixed[j],
                        heter_c_total[j],
                        lam1,
                        lambda2,
                        prior,
                        params.delta,
                    );
                    (j, kappa, beta, posterior, Some(lam1), Some(prior))
                }
            })
            .collect()
    };

    let mut hetersites: Vec<HeterSite> = Vec::with_capacity(heter_idx.len());
    for (j, kappa, beta, posterior, lam1, prior) in heter_estimates {
        let i = heter_idx[j];
        hetersites.push(HeterSite {
            pos: allsites.pos[i].clone(),
            motif: allsites.motif[i].clone(),
            motif_type: allsites.motif_type[i].clone(),
            treatment_x_rate: allsites.t_x_rate[i],
            control_x_rate: allsites.c_x_rate[i],
            treatment_fixed_rate: allsites.t_fixed_rate[i],
            control_fixed_rate: allsites.c_fixed_rate[i],
            treatment_fixed_count: allsites.t_fixed[i],
            control_fixed_count: allsites.c_fixed[i],
            treatment_depth: allsites.t_total[i],
            control_depth: allsites.c_total[i],
            beta_est: beta,
            kappa_est: kappa,
            binom_p: f64::NAN,
            binom_fdr: f64::NAN,
            fisher_p: f64::NAN,
            fisher_fdr: f64::NAN,
            log_lik_ratio: f64::NAN,
            lrt_p: f64::NAN,
            lrt_fdr: f64::NAN,
            lambda1: lam1,
            prior_methylated: prior,
            posterior: Some(posterior),
        });
    }

    Ok(PriorResult {
        homosites,
        hetersites,
        lambda1_summary: summary,
        lambda1_per_motif,
        lambda2,
    })
}
