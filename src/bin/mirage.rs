//! MIRAGE CLI: subcommands for compute-prior, estimate-empirical,
//! estimate-prior. All inputs/outputs are tab-separated tables.

use std::path::PathBuf;

use anyhow::{anyhow, Result};
use clap::{Parser, Subcommand, ValueEnum};

use mirage_rs::inference::{
    estimate_inference_with_empirical, estimate_inference_with_prior, BgMethod, BgTarget,
    EmpiricalParams, Lambda1Mode, NmerMode, PriorParams,
};
use mirage_rs::io::{
    read_count_table, read_motif_freq_bg, read_motif_freq_exp, write_tsv,
};
use mirage_rs::prior::compute_prior;

#[derive(Parser, Debug)]
#[command(name = "mirage", version, about = "MIRAGE Rust reimplementation")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Compute motif prior probabilities.
    ComputePrior {
        #[arg(long)]
        motif_freq_exp: PathBuf,
        #[arg(long)]
        motif_freq_bg: PathBuf,
        #[arg(long, default_value = "DRACH")]
        target_motif: String,
        #[arg(long, default_value_t = 3)]
        target_base_index: usize,
        #[arg(long, default_value_t = 0.007)]
        prob_methylation: f64,
        #[arg(long)]
        out: PathBuf,
    },
    /// Estimate site signal with empirical (constant) on-target rate.
    EstimateEmpirical {
        #[arg(long)]
        counts: PathBuf,
        #[arg(long, default_value_t = 0.001)]
        delta: f64,
        #[arg(long, default_value_t = 10)]
        depth_cutoff: u64,
        #[arg(long, default_value_t = 0.99)]
        homo_cutoff: f64,
        /// "auto", a numeric value, or path to a one-column file of site IDs.
        #[arg(long, default_value = "auto")]
        lambda1: String,
        #[arg(long, value_enum, default_value_t = CliBgMethod::Fisher)]
        bg_method: CliBgMethod,
        #[arg(long, value_enum, default_value_t = CliBgTarget::Treatment)]
        bg_target: CliBgTarget,
        #[arg(long, default_value_t = 0.95)]
        highly_methyl_cutoff: f64,
        #[arg(long)]
        seed: Option<u64>,
        #[arg(long, default_value_t = 1)]
        thread: usize,
        #[arg(long)]
        out_homo: PathBuf,
        #[arg(long)]
        out_heter: PathBuf,
        #[arg(long)]
        out_summary: PathBuf,
    },
    /// Estimate site signal with motif-aware prior.
    EstimatePrior {
        #[arg(long)]
        counts: PathBuf,
        #[arg(long)]
        prior: PathBuf,
        #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
        motif_specific: bool,
        #[arg(long, default_value_t = 0.001)]
        delta: f64,
        #[arg(long, default_value_t = 10)]
        depth_cutoff: u64,
        #[arg(long, default_value_t = 0.99)]
        homo_cutoff: f64,
        #[arg(long, value_enum, default_value_t = CliBgMethod::Fisher)]
        bg_method: CliBgMethod,
        #[arg(long, default_value_t = 0.95)]
        highly_methyl_cutoff: f64,
        #[arg(long)]
        seed: Option<u64>,
        #[arg(long, default_value_t = 1)]
        thread: usize,
        #[arg(long, default_value = "DRACH")]
        motif: String,
        #[arg(long, value_enum, default_value_t = CliNmer::Mer5)]
        nmer: CliNmer,
        #[arg(long)]
        out_homo: PathBuf,
        #[arg(long)]
        out_heter: PathBuf,
        #[arg(long)]
        out_summary: PathBuf,
    },
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum CliBgMethod {
    Binomial,
    Fisher,
    Lrt,
}

impl From<CliBgMethod> for BgMethod {
    fn from(v: CliBgMethod) -> Self {
        match v {
            CliBgMethod::Binomial => BgMethod::Binomial,
            CliBgMethod::Fisher => BgMethod::Fisher,
            CliBgMethod::Lrt => BgMethod::Lrt,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum CliBgTarget {
    Treatment,
    Control,
    Both,
}

impl From<CliBgTarget> for BgTarget {
    fn from(v: CliBgTarget) -> Self {
        match v {
            CliBgTarget::Treatment => BgTarget::Treatment,
            CliBgTarget::Control => BgTarget::Control,
            CliBgTarget::Both => BgTarget::Both,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum CliNmer {
    Mer5,
    F4,
    L4,
}

impl From<CliNmer> for NmerMode {
    fn from(v: CliNmer) -> Self {
        match v {
            CliNmer::Mer5 => NmerMode::Mer5,
            CliNmer::F4 => NmerMode::F4,
            CliNmer::L4 => NmerMode::L4,
        }
    }
}

fn parse_lambda1(s: &str) -> Result<Lambda1Mode> {
    if s == "auto" {
        return Ok(Lambda1Mode::Auto);
    }
    if let Ok(v) = s.parse::<f64>() {
        if (0.0..=1.0).contains(&v) {
            return Ok(Lambda1Mode::Fixed(v));
        }
    }
    // Otherwise treat as a path to a file with one site ID per line.
    let p = std::path::Path::new(s);
    if !p.exists() {
        return Err(anyhow!(
            "lambda1 must be 'auto', a numeric in [0,1], or path to site-ID list"
        ));
    }
    let content = std::fs::read_to_string(p)?;
    let ids: Vec<String> = content
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();
    Ok(Lambda1Mode::Site(ids))
}

#[derive(serde::Serialize)]
struct EmpSummaryRow {
    name: &'static str,
    value: f64,
}

#[derive(serde::Serialize)]
struct PriorSummaryRow {
    motif: String,
    lambda1: f64,
    tx_freq: f64,
    genome_freq: f64,
    freq: f64,
    motif_type: String,
    prior_methylated: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::ComputePrior {
            motif_freq_exp,
            motif_freq_bg,
            target_motif,
            target_base_index,
            prob_methylation,
            out,
        } => {
            let exp = read_motif_freq_exp(&motif_freq_exp)?;
            let bg = read_motif_freq_bg(&motif_freq_bg)?;
            let rows = compute_prior(
                &exp,
                &bg,
                &target_motif,
                target_base_index,
                prob_methylation,
            )?;
            write_tsv(&out, &rows)?;
        }
        Cmd::EstimateEmpirical {
            counts,
            delta,
            depth_cutoff,
            homo_cutoff,
            lambda1,
            bg_method,
            bg_target,
            highly_methyl_cutoff,
            seed,
            thread,
            out_homo,
            out_heter,
            out_summary,
        } => {
            if thread > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(thread)
                    .build_global()
                    .ok();
            }
            let recs = read_count_table(&counts)?;
            let params = EmpiricalParams {
                delta,
                depth_cutoff,
                homo_cutoff,
                lambda1: parse_lambda1(&lambda1)?,
                bg_method: bg_method.into(),
                bg_target: bg_target.into(),
                highly_methyl_cutoff,
                seed,
                thread,
            };
            let res = estimate_inference_with_empirical(&recs, &params)?;
            write_tsv(&out_homo, &res.homosites)?;
            write_tsv(&out_heter, &res.hetersites)?;
            let summary = vec![
                EmpSummaryRow {
                    name: "lambda1",
                    value: res.lambda1,
                },
                EmpSummaryRow {
                    name: "lambda2",
                    value: res.lambda2,
                },
            ];
            write_tsv(&out_summary, &summary)?;
        }
        Cmd::EstimatePrior {
            counts,
            prior,
            motif_specific,
            delta,
            depth_cutoff,
            homo_cutoff,
            bg_method,
            highly_methyl_cutoff,
            seed,
            thread,
            motif,
            nmer,
            out_homo,
            out_heter,
            out_summary,
        } => {
            if thread > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(thread)
                    .build_global()
                    .ok();
            }
            let recs = read_count_table(&counts)?;
            // Read prior table (already computed by compute_prior).
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_path(&prior)?;
            let mut prior_rows = Vec::new();
            for row in rdr.deserialize::<mirage_rs::prior::MotifPriorRow>() {
                prior_rows.push(row?);
            }
            let params = PriorParams {
                motif_specific,
                delta,
                depth_cutoff,
                homo_cutoff,
                bg_method: bg_method.into(),
                highly_methyl_cutoff,
                seed,
                thread,
                motif,
                nmer: nmer.into(),
            };
            let res = estimate_inference_with_prior(&recs, &params, &prior_rows)?;
            write_tsv(&out_homo, &res.homosites)?;
            write_tsv(&out_heter, &res.hetersites)?;
            let summary: Vec<PriorSummaryRow> = res
                .lambda1_summary
                .iter()
                .map(|r| PriorSummaryRow {
                    motif: r.motif.clone(),
                    lambda1: res
                        .lambda1_per_motif
                        .get(&r.motif)
                        .copied()
                        .unwrap_or(f64::NAN),
                    tx_freq: r.tx_freq,
                    genome_freq: r.genome_freq,
                    freq: r.freq,
                    motif_type: r.motif_type.clone(),
                    prior_methylated: r.prior_methylated,
                })
                .collect();
            write_tsv(&out_summary, &summary)?;
            eprintln!("lambda2 = {}", res.lambda2);
        }
    }
    Ok(())
}

