//! Motif prior probability table — port of `compute_prior.R`.

use std::collections::{BTreeSet, HashMap};

use anyhow::{anyhow, Result};

use crate::motif::{generate_motif_sequences, normalize_seq};

/// Per-motif observation frequency in the experimental signal set.
#[derive(Debug, Clone)]
pub struct MotifFreqExp {
    pub motif: String,
    pub freq: f64,
}

/// Per-motif counts in the transcriptome and genome.
#[derive(Debug, Clone)]
pub struct MotifFreqBg {
    pub motif: String,
    pub tx_count: f64,
    pub genome_count: f64,
}

/// Output row of `compute_prior`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MotifPriorRow {
    pub motif: String,
    pub tx_freq: f64,
    pub genome_freq: f64,
    pub freq: f64,
    pub motif_type: String,
    pub prior_methylated: f64,
}

fn motif_filter(s: &str, target_base: char, target_base_index: usize, motif_len: usize) -> bool {
    if s.len() != motif_len {
        return false;
    }
    let bytes = s.as_bytes();
    for (i, &c) in bytes.iter().enumerate() {
        let ch = c as char;
        if i + 1 == target_base_index {
            if ch != target_base {
                return false;
            }
        } else if !matches!(ch, 'A' | 'U' | 'C' | 'G') {
            return false;
        }
    }
    true
}

/// Direct port of `compute_prior`.
///
/// Behavior:
/// 1. Normalize motifs and target_motif: T -> U.
/// 2. Discard rows whose motif does not match `target_motif`'s structural
///    constraint (positions other than target_base_index must be ACGU; the
///    target position must equal `target_base`).
/// 3. Compute motif frequency tables `tx_freq = tx_count / sum(tx_count)`,
///    `genome_freq = genome_count / sum(genome_count)` over filtered bg rows.
/// 4. Tag each exp row by whether its motif is one of the concrete sequences
///    expanded from `target_motif`.
/// 5. Outer-merge bg and exp tables on `motif`.
/// 6. `prior_methylated = freq * prob_signal / tx_freq`.
pub fn compute_prior(
    motif_freq_exp: &[MotifFreqExp],
    motif_freq_bg: &[MotifFreqBg],
    target_motif: &str,
    target_base_index: usize,
    prob_signal: f64,
) -> Result<Vec<MotifPriorRow>> {
    let target_motif = normalize_seq(target_motif);
    let motif_len = target_motif.chars().count();
    if target_base_index == 0 || target_base_index > motif_len {
        return Err(anyhow!(
            "target_base_index ({}) must be between 1 and the motif length ({})",
            target_base_index,
            motif_len
        ));
    }
    let target_base = target_motif
        .chars()
        .nth(target_base_index - 1)
        .ok_or_else(|| anyhow!("target_base_index out of range"))?;
    if !matches!(target_base, 'A' | 'C' | 'U' | 'G') {
        return Err(anyhow!("The target base must be one of A, C, U or G!"));
    }

    let exp: Vec<MotifFreqExp> = motif_freq_exp
        .iter()
        .filter_map(|r| {
            let m = normalize_seq(&r.motif);
            if motif_filter(&m, target_base, target_base_index, motif_len) {
                Some(MotifFreqExp { motif: m, freq: r.freq })
            } else {
                None
            }
        })
        .collect();

    let bg: Vec<MotifFreqBg> = motif_freq_bg
        .iter()
        .filter_map(|r| {
            let m = normalize_seq(&r.motif);
            if motif_filter(&m, target_base, target_base_index, motif_len) {
                Some(MotifFreqBg {
                    motif: m,
                    tx_count: r.tx_count,
                    genome_count: r.genome_count,
                })
            } else {
                None
            }
        })
        .collect();

    let sum_tx: f64 = bg.iter().map(|r| r.tx_count).sum();
    let sum_genome: f64 = bg.iter().map(|r| r.genome_count).sum();

    let bg_index: HashMap<String, &MotifFreqBg> =
        bg.iter().map(|r| (r.motif.clone(), r)).collect();
    let exp_index: HashMap<String, &MotifFreqExp> =
        exp.iter().map(|r| (r.motif.clone(), r)).collect();

    let target_seqs: BTreeSet<String> = generate_motif_sequences(&target_motif)?
        .into_iter()
        .collect();

    let mut all_motifs: BTreeSet<String> = BTreeSet::new();
    for k in bg_index.keys() {
        all_motifs.insert(k.clone());
    }
    for k in exp_index.keys() {
        all_motifs.insert(k.clone());
    }

    let non_target_label = format!("non{}", target_motif);
    let mut out: Vec<MotifPriorRow> = Vec::with_capacity(all_motifs.len());
    for m in all_motifs {
        let bg_row = bg_index.get(&m);
        let exp_row = exp_index.get(&m);
        let tx_freq = bg_row
            .map(|r| if sum_tx > 0.0 { r.tx_count / sum_tx } else { f64::NAN })
            .unwrap_or(f64::NAN);
        let genome_freq = bg_row
            .map(|r| {
                if sum_genome > 0.0 {
                    r.genome_count / sum_genome
                } else {
                    f64::NAN
                }
            })
            .unwrap_or(f64::NAN);
        let freq = exp_row.map(|r| r.freq).unwrap_or(f64::NAN);
        let motif_type = if exp_row.is_some() {
            if target_seqs.contains(&m) {
                target_motif.clone()
            } else {
                non_target_label.clone()
            }
        } else {
            String::new()
        };
        let prior_methylated = if tx_freq.is_nan() || tx_freq == 0.0 || freq.is_nan() {
            f64::NAN
        } else {
            freq * prob_signal / tx_freq
        };
        out.push(MotifPriorRow {
            motif: m,
            tx_freq,
            genome_freq,
            freq,
            motif_type,
            prior_methylated,
        });
    }
    Ok(out)
}
