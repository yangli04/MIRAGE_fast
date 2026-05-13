#!/usr/bin/env Rscript
# Compare R and Rust outputs side-by-side. Reports max/median absolute
# differences in numeric columns and counts of mismatched site IDs.

R_dir <- "validation/R_out"
Rust_dir <- "validation/rust_out"

read_tsv <- function(p) read.delim(p, stringsAsFactors = FALSE)

approx_equal <- function(a, b, tol = 1e-6) {
  na_a <- is.na(a); na_b <- is.na(b)
  if (any(na_a != na_b)) return(FALSE)
  ok <- !na_a
  all(abs(a[ok] - b[ok]) < tol | abs((a[ok] - b[ok]) /
        pmax(abs(b[ok]), 1e-12)) < tol)
}

summarize_diff <- function(a, b, label) {
  if (length(a) != length(b)) {
    cat(sprintf("[%s] LENGTH MISMATCH: R=%d, Rust=%d\n",
                label, length(a), length(b)))
    return(invisible(NULL))
  }
  na_a <- is.na(a); na_b <- is.na(b)
  if (any(na_a != na_b)) {
    cat(sprintf("[%s] NA pattern mismatch: R has %d, Rust has %d\n",
                label, sum(na_a), sum(na_b)))
  }
  ok <- !na_a & !na_b
  if (!any(ok)) {
    cat(sprintf("[%s] all NA in one side\n", label))
    return(invisible(NULL))
  }
  d <- abs(a[ok] - b[ok])
  cat(sprintf("[%s] n=%d max|diff|=%.3e median|diff|=%.3e\n",
              label, sum(ok), max(d), median(d)))
}

cat("=== compute_prior ===\n")
prior_R <- read_tsv(file.path(R_dir, "prior.tsv"))
prior_rust <- read_tsv(file.path(Rust_dir, "prior.tsv"))
prior_R <- prior_R[order(prior_R$motif), ]
prior_rust <- prior_rust[order(prior_rust$motif), ]
cat("R rows:", nrow(prior_R), " Rust rows:", nrow(prior_rust), "\n")
common <- intersect(prior_R$motif, prior_rust$motif)
ar <- prior_R[match(common, prior_R$motif), ]
br <- prior_rust[match(common, prior_rust$motif), ]
summarize_diff(ar$tx_freq, br$tx_freq, "prior$tx_freq")
summarize_diff(ar$genome_freq, br$genome_freq, "prior$genome_freq")
summarize_diff(ar$freq, br$freq, "prior$freq")
summarize_diff(ar$prior_methylated, br$prior_methylated, "prior$prior_methylated")

cat("\n=== empirical: lambdas ===\n")
sumR <- readLines(file.path(R_dir, "emp_summary.tsv"))
sumRust <- read_tsv(file.path(Rust_dir, "emp_summary.tsv"))
cat("R:\n"); print(sumR)
cat("Rust:\n"); print(sumRust)

cat("\n=== empirical: homo sites ===\n")
homo_R <- read_tsv(file.path(R_dir, "emp_homo.tsv"))
homo_rust <- read_tsv(file.path(Rust_dir, "emp_homo.tsv"))
cat("R rows:", nrow(homo_R), " Rust rows:", nrow(homo_rust), "\n")
homo_R <- homo_R[order(homo_R$pos), ]
homo_rust <- homo_rust[order(homo_rust$pos), ]
stopifnot(identical(homo_R$pos, homo_rust$pos))
summarize_diff(homo_R$beta_est, homo_rust$beta_est, "homo$beta_est")
summarize_diff(homo_R$binom_p, homo_rust$binom_p, "homo$binom_p")
summarize_diff(homo_R$fisher_p, homo_rust$fisher_p, "homo$fisher_p")
summarize_diff(homo_R$lrt_p, homo_rust$lrt_p, "homo$lrt_p")
summarize_diff(homo_R$lrt_fdr, homo_rust$lrt_fdr, "homo$lrt_fdr")

cat("\n=== empirical: heter sites ===\n")
het_R <- read_tsv(file.path(R_dir, "emp_heter.tsv"))
het_rust <- read_tsv(file.path(Rust_dir, "emp_heter.tsv"))
cat("R rows:", nrow(het_R), " Rust rows:", nrow(het_rust), "\n")
het_R <- het_R[order(het_R$pos), ]
het_rust <- het_rust[order(het_rust$pos), ]
stopifnot(identical(het_R$pos, het_rust$pos))
summarize_diff(het_R$beta_est, het_rust$beta_est, "heter$beta_est")
summarize_diff(het_R$kappa_est, het_rust$kappa_est, "heter$kappa_est")
summarize_diff(het_R$binom_p, het_rust$binom_p, "heter$binom_p")
summarize_diff(het_R$fisher_p, het_rust$fisher_p, "heter$fisher_p")
summarize_diff(het_R$lrt_p, het_rust$lrt_p, "heter$lrt_p")
summarize_diff(het_R$lrt_fdr, het_rust$lrt_fdr, "heter$lrt_fdr")

cat("\n=== prior: homo posterior ===\n")
phomo_R <- read_tsv(file.path(R_dir, "prior_homo.tsv"))
phomo_rust <- read_tsv(file.path(Rust_dir, "prior_homo.tsv"))
cat("R rows:", nrow(phomo_R), " Rust rows:", nrow(phomo_rust), "\n")
phomo_R <- phomo_R[order(phomo_R$pos), ]
phomo_rust <- phomo_rust[order(phomo_rust$pos), ]
common <- intersect(phomo_R$pos, phomo_rust$pos)
ar <- phomo_R[match(common, phomo_R$pos), ]
br <- phomo_rust[match(common, phomo_rust$pos), ]
summarize_diff(ar$beta_est, br$beta_est, "prior_homo$beta_est")
summarize_diff(ar$posterior, br$posterior, "prior_homo$posterior")

cat("\n=== prior: heter posterior ===\n")
phet_R <- read_tsv(file.path(R_dir, "prior_heter.tsv"))
phet_rust <- read_tsv(file.path(Rust_dir, "prior_heter.tsv"))
cat("R rows:", nrow(phet_R), " Rust rows:", nrow(phet_rust), "\n")
phet_R <- phet_R[order(phet_R$pos), ]
phet_rust <- phet_rust[order(phet_rust$pos), ]
common <- intersect(phet_R$pos, phet_rust$pos)
ar <- phet_R[match(common, phet_R$pos), ]
br <- phet_rust[match(common, phet_rust$pos), ]
summarize_diff(ar$beta_est, br$beta_est, "prior_heter$beta_est")
summarize_diff(ar$kappa_est, br$kappa_est, "prior_heter$kappa_est")
summarize_diff(ar$posterior, br$posterior, "prior_heter$posterior")

cat("\n=== Done. ===\n")
