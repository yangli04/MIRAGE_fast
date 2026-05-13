//! TSV I/O helpers for MIRAGE.

use std::path::Path;

use anyhow::{anyhow, Context, Result};
use serde::Deserialize;

use crate::inference::SiteRecord;
use crate::prior::{MotifFreqBg, MotifFreqExp};

#[derive(Debug, Deserialize)]
struct CountRow {
    pos: String,
    motif: String,
    #[serde(rename = "type")]
    motif_type: String,
    treated_fixed_count: u64,
    treated_depth: u64,
    control_fixed_count: u64,
    control_depth: u64,
}

/// Read a count table with columns:
/// pos, motif, type, treated_fixed_count, treated_depth, control_fixed_count, control_depth
pub fn read_count_table(path: &Path) -> Result<Vec<SiteRecord>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let mut out = Vec::new();
    for row in rdr.deserialize::<CountRow>() {
        let r = row?;
        out.push(SiteRecord {
            pos: r.pos,
            motif: r.motif,
            motif_type: r.motif_type,
            treated_fixed_count: r.treated_fixed_count,
            treated_depth: r.treated_depth,
            control_fixed_count: r.control_fixed_count,
            control_depth: r.control_depth,
        });
    }
    if out.is_empty() {
        return Err(anyhow!("count table is empty: {}", path.display()));
    }
    Ok(out)
}

#[derive(Debug, Deserialize)]
struct ExpRow {
    motif: String,
    freq: f64,
}

#[derive(Debug, Deserialize)]
struct BgRow {
    motif: String,
    tx_count: f64,
    genome_count: f64,
}

pub fn read_motif_freq_exp(path: &Path) -> Result<Vec<MotifFreqExp>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let mut out = Vec::new();
    for row in rdr.deserialize::<ExpRow>() {
        let r = row?;
        out.push(MotifFreqExp {
            motif: r.motif,
            freq: r.freq,
        });
    }
    Ok(out)
}

pub fn read_motif_freq_bg(path: &Path) -> Result<Vec<MotifFreqBg>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let mut out = Vec::new();
    for row in rdr.deserialize::<BgRow>() {
        let r = row?;
        out.push(MotifFreqBg {
            motif: r.motif,
            tx_count: r.tx_count,
            genome_count: r.genome_count,
        });
    }
    Ok(out)
}

/// Write a slice of serializable rows out as a tab-separated file.
pub fn write_tsv<T: serde::Serialize>(path: &Path, rows: &[T]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("creating {}", path.display()))?;
    for r in rows {
        wtr.serialize(r)?;
    }
    wtr.flush()?;
    Ok(())
}
