//! Degenerate base expansion and motif sequence generation.

use anyhow::{anyhow, Result};

/// Map a degenerate IUPAC base symbol to its expanded RNA base set.
/// T is normalized to U upstream.
pub fn degenerate_bases(b: char) -> Result<&'static [char]> {
    Ok(match b {
        'A' => &['A'],
        'C' => &['C'],
        'G' => &['G'],
        'U' => &['U'],
        'R' => &['A', 'G'],
        'Y' => &['C', 'U'],
        'S' => &['G', 'C'],
        'W' => &['A', 'U'],
        'K' => &['G', 'U'],
        'M' => &['A', 'C'],
        'B' => &['C', 'G', 'U'],
        'D' => &['A', 'G', 'U'],
        'H' => &['A', 'C', 'U'],
        'V' => &['A', 'C', 'G'],
        'N' => &['A', 'C', 'G', 'U'],
        other => return Err(anyhow!("unknown degenerate base: {}", other)),
    })
}

/// Normalize a motif/sequence by uppercasing and converting T to U.
pub fn normalize_seq(s: &str) -> String {
    s.chars()
        .map(|c| match c.to_ascii_uppercase() {
            'T' => 'U',
            x => x,
        })
        .collect()
}

/// Expand a degenerate motif (e.g. "DRACH") into all concrete sequences.
/// Returns sequences as strings in the same lexicographic order as
/// `expand.grid` in R (first column varies fastest).
pub fn generate_motif_sequences(motif: &str) -> Result<Vec<String>> {
    let motif = normalize_seq(motif);
    let bases: Vec<&'static [char]> = motif
        .chars()
        .map(degenerate_bases)
        .collect::<Result<_>>()?;

    // Match R `expand.grid` ordering: first column varies fastest.
    // We build sequences position-by-position; at each position cycling
    // through that position's bases. Using "first varies fastest" means
    // for position i, the cycle length before repeating is product of
    // sizes of positions 0..i. We replicate this by iterative product
    // updating column-wise.
    let n = bases.len();
    let total: usize = bases.iter().map(|b| b.len()).product();
    let mut cols: Vec<Vec<char>> = vec![Vec::with_capacity(total); n];
    let mut cycle = 1usize;
    for i in 0..n {
        let len_i = bases[i].len();
        let block = cycle * len_i;
        for row in 0..total {
            let idx = (row % block) / cycle;
            cols[i].push(bases[i][idx]);
        }
        cycle = block;
    }
    let seqs: Vec<String> = (0..total)
        .map(|r| (0..n).map(|c| cols[c][r]).collect())
        .collect();
    Ok(seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nnunn_has_256() {
        let seqs = generate_motif_sequences("NNUNN").unwrap();
        assert_eq!(seqs.len(), 256);
        for s in &seqs {
            assert_eq!(s.chars().nth(2), Some('U'));
        }
    }

    #[test]
    fn drach_has_18() {
        let seqs = generate_motif_sequences("DRACH").unwrap();
        assert_eq!(seqs.len(), 3 * 2 * 1 * 1 * 3);
    }

    #[test]
    fn first_varies_fastest() {
        let seqs = generate_motif_sequences("RY").unwrap();
        // R = A,G; Y = C,U. expand.grid gives:
        // AC GC AU GU  (column 1 cycles fastest)
        assert_eq!(seqs, vec!["AC", "GC", "AU", "GU"]);
    }

    #[test]
    fn t_to_u() {
        assert_eq!(normalize_seq("ACGT"), "ACGU");
    }
}
