# mirage_rs

Rust port of [MIRAGE](https://github.com/yangli04/MIRAGE) — Mutation-encoded
Inference of RNA Activity via Generative Effects.

`mirage_rs` estimates site-level RNA mutation or conversion signal from matched
treatment/control read-count tables. It provides the same three top-level
analyses as the R package, with per-site numerical results that match the R
reference to within Brent tolerance (≤ 1e-7) on the canonical PAR-CLIP
example.

| R function                           | Rust API / CLI |
|---|---|
| `compute_prior()`                    | `mirage_rs::compute_prior` / `mirage compute-prior` |
| `estimate_inference_with_empirical()`| `mirage_rs::estimate_inference_with_empirical` / `mirage estimate-empirical` |
| `estimate_inference_with_prior()`    | `mirage_rs::estimate_inference_with_prior` / `mirage estimate-prior` |

---

## Build

Requires Rust ≥ 1.75 (any recent stable will work; tested on 1.93).

```bash
cargo build --release
```

This produces:

- `target/release/libmirage_rs.rlib` — Rust library
- `target/release/mirage` — CLI

Run the unit tests:

```bash
cargo test --release
```

---

## CLI quick start

The CLI reads tab-separated tables and writes tab-separated tables. All the
examples below assume you are running from the project root and that the
PAR-CLIP example has been copied to `validation/parclip.tsv` (see
`validation/run_R_reference.R` for one way to produce it).

### 1. Empirical inference (one global on-target rate)

```bash
./target/release/mirage estimate-empirical \
  --counts validation/parclip.tsv \
  --delta 0.001 \
  --depth-cutoff 10 \
  --homo-cutoff 0.99 \
  --lambda1 auto \
  --bg-method lrt \
  --bg-target both \
  --highly-methyl-cutoff 0.95 \
  --seed 123 \
  --thread 4 \
  --out-homo  out/emp_homo.tsv \
  --out-heter out/emp_heter.tsv \
  --out-summary out/emp_summary.tsv
```

`--lambda1` accepts:

- `auto` — estimate from the upper quantile of candidate signal sites
  (default).
- a number in `[0, 1]` — fix `lambda1` explicitly.
- a path to a one-column file of site IDs — estimate from that user-defined
  high-signal set.

### 2. Compute motif priors

```bash
./target/release/mirage compute-prior \
  --motif-freq-exp motif_freq_exp.tsv \
  --motif-freq-bg  motif_freq_bg.tsv \
  --target-motif NNUNN \
  --target-base-index 3 \
  --prob-methylation 0.01 \
  --out out/prior.tsv
```

### 3. Prior-aware inference

```bash
./target/release/mirage estimate-prior \
  --counts validation/parclip.tsv \
  --prior  out/prior.tsv \
  --motif-specific true \
  --motif NNUNN \
  --nmer mer5 \
  --delta 0.001 \
  --depth-cutoff 10 \
  --bg-method lrt \
  --highly-methyl-cutoff 0.95 \
  --seed 123 \
  --thread 4 \
  --out-homo  out/prior_homo.tsv \
  --out-heter out/prior_heter.tsv \
  --out-summary out/prior_summary.tsv
```

`--nmer` is one of:

- `mer5` — estimate `lambda1` separately for every concrete 5-mer.
- `f4`   — group concrete 5-mers by their first 4 bases.
- `l4`   — group concrete 5-mers by their last 4 bases.

`--motif-specific false` falls back to the single overall `lambda1` for every
target motif.

`-h` / `--help` lists every flag for any subcommand.

---

## Input file formats

All inputs are TSV with a header row.

### Count table — `--counts`

| column                  | type   | description |
|---|---|---|
| `pos`                   | string | site ID, e.g. `chr1_14404_-` |
| `motif`                 | string | local motif sequence (e.g. `AAACC`); `T` is treated as `U` |
| `type`                  | string | motif class (free-form, e.g. `NNUNN`, `DRACH`, `nonDRACH`) |
| `treated_fixed_count`   | u64    | non-mutated read count, treatment |
| `treated_depth`         | u64    | total read depth, treatment |
| `control_fixed_count`   | u64    | non-mutated read count, control |
| `control_depth`         | u64    | total read depth, control |

### Motif frequency tables — `compute-prior`

`--motif-freq-exp`:

| column | type   | description |
|---|---|---|
| `motif`| string | 5-mer motif sequence |
| `freq` | f64    | motif frequency among independently defined signal sites |

`--motif-freq-bg`:

| column         | type   | description |
|---|---|---|
| `motif`        | string | 5-mer motif sequence |
| `tx_count`     | f64    | transcriptome motif count |
| `genome_count` | f64    | genome motif count |

### Prior table — `estimate-prior --prior`

The TSV produced by `mirage compute-prior`, with columns `motif`, `tx_freq`,
`genome_freq`, `freq`, `motif_type`, `prior_methylated`. You can also build
this externally as long as it has those columns.

---

## Output file formats

### `out_homo` / `out_heter`

One row per site with the input fields plus per-site estimates. Columns:

| column | description |
|---|---|
| `pos`, `motif`, `motif_type` | from input |
| `treatment_x_rate`, `control_x_rate` | observed mutation rate (`1 − fixed/depth`) |
| `treatment_fixed_rate`, `control_fixed_rate` | observed non-mutation rate |
| `treatment_fixed_count`, `control_fixed_count`, `treatment_depth`, `control_depth` | input counts |
| `binom_p`, `binom_fdr` | one-sided binomial test against the inferred background rate (BH-adjusted) |
| `fisher_p`, `fisher_fdr` | one-sided Fisher's exact test (BH-adjusted) |
| `log_lik_ratio`, `lrt_p`, `lrt_fdr` | methylSig-style LRT (BH-adjusted) |
| `beta_est` | inferred site-level signal |
| `kappa_est` (heter only) | allele-fraction estimate |
| `lambda1`, `prior_methylated`, `posterior` | populated by `estimate-prior` only |

Statistical-test columns are written as `NaN` by `estimate-prior` (R behavior:
priors-mode does not run a per-site test step).

### `out_summary`

For `estimate-empirical`: two rows, `lambda1` and `lambda2`.

For `estimate-prior`: per-motif summary with columns `motif`, `lambda1`,
`tx_freq`, `genome_freq`, `freq`, `motif_type`, `prior_methylated`. The
overall `lambda2` is logged to stderr.

---

## Library usage (Rust)

Add a path dependency and call the high-level functions directly.

```toml
[dependencies]
mirage_rs = { path = "/data/yangli/20260509_rust_implement_methylhunter" }
```

### Empirical inference

```rust
use mirage_rs::{
    estimate_inference_with_empirical, BgMethod, BgTarget,
    EmpiricalParams, Lambda1Mode, SiteRecord,
};

let records = vec![
    SiteRecord {
        pos: "chr1_14404_-".into(),
        motif: "AAACC".into(),
        motif_type: "NNUNN".into(),
        treated_fixed_count: 16,
        treated_depth: 22,
        control_fixed_count: 22,
        control_depth: 22,
    },
    // ... more sites ...
];

let params = EmpiricalParams {
    delta: 0.001,
    depth_cutoff: 10,
    homo_cutoff: 0.99,
    lambda1: Lambda1Mode::Auto,
    bg_method: BgMethod::Lrt,
    bg_target: BgTarget::Both,
    highly_methyl_cutoff: 0.95,
    seed: Some(123),
    thread: 4,
};

let res = estimate_inference_with_empirical(&records, &params)?;
println!("lambda1 = {}, lambda2 = {}", res.lambda1, res.lambda2);
for site in &res.homosites {
    println!("{}\t{}\t{}", site.pos, site.beta_est, site.lrt_fdr);
}
```

### Compute prior + prior-aware inference

```rust
use mirage_rs::{
    compute_prior, estimate_inference_with_prior,
    BgMethod, MotifFreqBg, MotifFreqExp, NmerMode, PriorParams,
};

let exp = vec![MotifFreqExp { motif: "AAACC".into(), freq: 0.01 }, /* ... */];
let bg  = vec![MotifFreqBg  { motif: "AAACC".into(), tx_count: 100.0, genome_count: 120.0 }, /* ... */];

let prior = compute_prior(&exp, &bg, "DRACH", 3, 0.007)?;

let params = PriorParams {
    motif_specific: true,
    delta: 0.001,
    depth_cutoff: 10,
    homo_cutoff: 0.99,
    bg_method: BgMethod::Lrt,
    highly_methyl_cutoff: 0.95,
    seed: Some(123),
    thread: 4,
    motif: "DRACH".into(),
    nmer: NmerMode::Mer5,
};

let res = estimate_inference_with_prior(&records, &params, &prior)?;
println!("lambda2 = {}", res.lambda2);
for site in &res.homosites {
    println!("{}\t{:?}\t{:?}", site.pos, site.beta_est, site.posterior);
}
```

### Lower-level helpers

If you only need a sub-component, the underlying modules are public:

| Module | Public items |
|---|---|
| `mirage_rs::motif`     | `generate_motif_sequences`, `degenerate_bases`, `normalize_seq` |
| `mirage_rs::stats`     | `binom_test_greater`, `fisher_test_greater`, `log_lik_ratio_test`, `bh_adjust`, `binom_pmf`, `binom_cdf`, `binom_sf_ge` |
| `mirage_rs::optim`     | `brent_min`, `quantile_type7`, `median`, `nelder_mead_2d` |
| `mirage_rs::mle`       | `mle_for_lambda1`, `mle_for_lambda2`, `mle_for_beta_homo`, `mle_for_beta_heter`, `mle_for_kappa_heter`, `mle_joint_beta_kappa_heter`, `bayesian_inference_beta_homo`, `bayesian_inference_beta_kappa_heter` |
| `mirage_rs::prior`     | `compute_prior`, `MotifFreqExp`, `MotifFreqBg`, `MotifPriorRow` |
| `mirage_rs::inference` | `EmpiricalParams`, `PriorParams`, `EmpiricalResult`, `PriorResult`, `SiteRecord`, `HomoSite`, `HeterSite`, plus the two top-level functions |
| `mirage_rs::io`        | `read_count_table`, `read_motif_freq_exp`, `read_motif_freq_bg`, `write_tsv` |

---

## Parameters at a glance

The defaults match the R reference one-for-one.

| Parameter | Default | Meaning |
|---|---|---|
| `delta` | 0.001 | sequencing error rate |
| `depth_cutoff` | 10 | minimum treatment & control depth to keep a site |
| `homo_cutoff` | 0.99 | control non-mutation-rate threshold for "homozygous" |
| `bg_method` | `fisher` | which test selects background sites (`binomial`, `fisher`, `lrt`) |
| `bg_target` | `treatment` | which sample's counts feed `lambda2` MLE (`treatment`, `control`, `both`) |
| `highly_methyl_cutoff` | 0.95 | upper-quantile cutoff for the high-signal set |
| `lambda1` | `auto` | `auto` / numeric / site list |
| `seed` | none | RNG seed for the lambda2 subsample |
| `thread` | 1 | rayon worker threads |
| `motif` | `DRACH` | target motif (degenerate IUPAC; `T` ↔ `U`) |
| `nmer` | `mer5` | motif grouping for per-motif `lambda1` (`mer5`, `f4`, `l4`) |
| `motif_specific` | `true` | per-motif vs. single overall `lambda1` |
| `target_base_index` | 3 | 1-based position of the target base inside `target_motif` |
| `prob_methylation` | 0.007 | overall prior probability of true signal at the target base |

Both the CLI flags (kebab-case) and the Rust struct fields (snake_case) use
these names.

---

## Validation against R MIRAGE

`validation/run_R_reference.R` runs the R reference on the bundled
`pos.read.count.example.parclip.txt` and writes its outputs to
`validation/R_out/`. `validation/compare.R` then re-reads the matching files
under `validation/rust_out/` and reports max/median absolute differences.

```bash
Rscript validation/run_R_reference.R    # writes validation/R_out/
# then run the three CLI commands writing into validation/rust_out/, e.g.
./target/release/mirage estimate-empirical --counts validation/parclip.tsv \
    --bg-method lrt --bg-target both --seed 123 \
    --out-homo validation/rust_out/emp_homo.tsv \
    --out-heter validation/rust_out/emp_heter.tsv \
    --out-summary validation/rust_out/emp_summary.tsv
Rscript validation/compare.R
```

Observed agreement on the PAR-CLIP example (989 homozygous + 10 heterozygous
sites):

| Quantity | Max abs diff Rust vs R |
|---|---|
| `compute_prior` (`tx_freq`, `genome_freq`, `prior_methylated`) | 1e-17 (machine ε) |
| `compute_prior$freq` | 0 |
| empirical `lambda1`, `lambda2` | identical |
| empirical homo `beta_est` | 2.5e-8 (Brent tolerance) |
| empirical homo p-values & FDRs | < 5e-8 |
| empirical heter `beta_est` / `kappa_est` | 1e-5 (R rounds to 5 decimals) |
| prior homo `beta_est` / `posterior` | < 5e-8 |
| prior heter `beta_est` / `kappa_est` / `posterior` | < 1e-8 |

Site-set sizes (homozygous, background, high-signal candidates,
heterozygous) match R exactly.

---

## Run-time benchmark

Same machine, single config: `bg_method = lrt`, `bg_target = both`,
`lambda1 = auto`, `seed = 123`. R timing uses `data.table::fread` and is
restricted to the same call you would have run from the R package; no
`compute_prior` is in the timed region. Memory is `max-rss` reported by
`/usr/bin/time`. Each run was repeated and the warm number is reported. The
synthetic 100k dataset is the bundled 999-site PAR-CLIP example tiled 100×
with `pos` disambiguated.

### `estimate_inference_with_empirical`

| Dataset       | Implementation     | Wall time | Max RSS  | Speed-up vs R |
|---|---|---:|---:|---:|
| 999 sites     | R reference, t=1   | 1.38 s    | 114 MB   | 1×            |
| 999 sites     | Rust, t=1          | 0.04 s    | 4 MB     | **35×**       |
| 999 sites     | Rust, t=4          | 0.03 s    | 4 MB     | 46×           |
| 999 sites     | Rust, t=8          | 0.03 s    | 4 MB     | 46×           |
| 99,900 sites  | R reference, t=1   | 71.3 s    | 211 MB   | 1×            |
| 99,900 sites  | Rust, t=1          | 1.50 s    | 98 MB    | **48×**       |
| 99,900 sites  | Rust, t=4          | 1.45 s    | 97 MB    | 49×           |
| 99,900 sites  | Rust, t=8          | 1.41 s    | 95 MB    | 51×           |

### `estimate_inference_with_prior` (motif_specific = FALSE)

| Dataset       | Implementation     | Wall time | Max RSS  | Speed-up vs R |
|---|---|---:|---:|---:|
| 99,900 sites  | R reference, t=1   | 50.9 s    | 232 MB   | 1×            |
| 99,900 sites  | Rust, t=1          | 1.40 s    | 107 MB   | **36×**       |
| 99,900 sites  | Rust, t=4          | 1.34 s    | 104 MB   | 38×           |
| 99,900 sites  | Rust, t=8          | 1.28 s    | 106 MB   | 40×           |

### What to read from the numbers

- **Single-thread Rust is ~35–50× faster than R** across both pipelines and
  both dataset sizes, with ~2× lower memory footprint on the larger run.
- **Threading saturates fast on the bundled examples.** On the PAR-CLIP
  dataset the Brent calls take a few microseconds each, so rayon's overhead
  dominates beyond `--thread 1`. Threading still helps for the per-site
  binomial/Fisher/LRT passes when those iterate over hundreds of thousands
  of sites — useful in the millions-of-sites regime closer to a full PAR-CLIP
  run.
- **R startup (~0.5 s) is amortized.** R's `Rscript` paid roughly 0.5 s of
  fixed startup, so for the 999-site case the inner R timer reports closer
  to 0.8 s. Even with that correction, single-thread Rust is ~25× faster on
  the small set and the gap widens linearly with site count.
- **Reproduce locally** with `bash validation/bench.sh`. The script builds
  release artifacts, generates the synthetic 100k dataset on first run, and
  prints `/usr/bin/time -v`-style measurements for every cell of the table.

---

## Notes / behavioral parity with R

- **T ↔ U.** All motif strings are uppercased and `T` is converted to `U`.
- **Quantile.** R's default `quantile()` (type 7, continuous interpolation) is
  used everywhere — see `optim::quantile_type7`.
- **BH FDR.** Implements R's `p.adjust(method = "BH")` exactly, including the
  cumulative-min step.
- **`fisher.test(alternative = "greater")`.** Implemented as the upper-tail
  hypergeometric sum; numerically equal to R's exact computation for the
  small marginals encountered in MIRAGE.
- **Joint `(beta, kappa)` MLE on heter sites.** R uses L-BFGS-B and rounds the
  result to 5 decimals. The Rust port runs alternating Brent on each
  coordinate (coordinate descent), which is more reliable than 2-D
  Nelder-Mead on this likelihood and reproduces R's rounded output.
- **Empty high-signal set.** When `MLE_for_lambda1` is called on an empty
  count vector, the Rust port returns the initial value (`0.2`), matching
  R's Brent behavior on a constant likelihood. Truly zero-signal datasets
  produce a meaningless `lambda1` either way; both implementations propagate
  it identically downstream.
- **`estimate_inference_with_prior` filters.** Rows are kept iff their motif
  has a non-NaN `prior_methylated`; rows whose motif has NaN per-motif
  `lambda1` (under-powered groups in motif-specific mode) are kept with NaN
  `beta_est` / `posterior`, matching the R output shape.

---

## License

MIT, matching the upstream R MIRAGE package.
