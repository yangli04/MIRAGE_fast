#!/usr/bin/env Rscript
# Run the R MIRAGE reference on the bundled example, write outputs as TSV
# for comparison with the Rust port.

suppressPackageStartupMessages({
  library(MIRAGE)
})

set.seed(123)

count_path <- system.file("extdata", "pos.read.count.example.parclip.txt",
                          package = "MIRAGE")
count_table <- read.delim(count_path, stringsAsFactors = FALSE)

out_dir <- file.path("validation", "R_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. compute_prior on the table itself (mirrors compute_prior example).
motif_counts <- table(count_table$motif)
motif_freq_exp <- data.frame(
  motif = names(motif_counts),
  freq = as.numeric(motif_counts) / sum(motif_counts),
  stringsAsFactors = FALSE
)
motif_freq_bg <- data.frame(
  motif = names(motif_counts),
  tx_count = as.numeric(motif_counts),
  genome_count = as.numeric(motif_counts),
  stringsAsFactors = FALSE
)
write.table(motif_freq_exp,
            file = file.path(out_dir, "motif_freq_exp.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(motif_freq_bg,
            file = file.path(out_dir, "motif_freq_bg.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
prior <- compute_prior(motif_freq_exp, motif_freq_bg,
                       target_motif = "NNUNN",
                       target_base_index = 3,
                       prob_methylation = 0.01)
write.table(prior, file = file.path(out_dir, "prior.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 2. estimate_inference_with_empirical
emp <- with(count_table, estimate_inference_with_empirical(
  cbind(pos, motif, type),
  treated_fixed_count, control_fixed_count,
  treated_depth, control_depth,
  delta = 0.001, depth.cutoff = 10,
  lambda1 = "auto", top.sites = NULL,
  bg.method = "lrt", bg.target = "both",
  highly.methyl.cutoff = 0.95,
  seed = 123, thread = 1
))

cat("R lambda1 =", emp$lambda1, "\n")
cat("R lambda2 =", emp$lambda2, "\n")

write.table(emp$homosites,
            file = file.path(out_dir, "emp_homo.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(emp$hetersites,
            file = file.path(out_dir, "emp_heter.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c(paste0("lambda1\t", emp$lambda1),
             paste0("lambda2\t", emp$lambda2)),
           con = file.path(out_dir, "emp_summary.tsv"))

# 3. estimate_inference_with_prior (motif_specific = FALSE for stable comparison)
pri <- with(count_table, estimate_inference_with_prior(
  cbind(pos, motif, type),
  treated_fixed_count, control_fixed_count,
  treated_depth, control_depth,
  motif_specific = FALSE,
  delta = 0.001, depth.cutoff = 10,
  bg.method = "lrt",
  highly.methyl.cutoff = 0.95,
  ref.freq.tab = prior,
  seed = 123, thread = 1,
  motif = "NNUNN", Nmer = "5mer"
))
cat("R prior lambda2 =", pri$lambda2, "\n")
cat("R prior lambda1 (head):\n")
print(head(pri$lambda1))

write.table(pri$homosites,
            file = file.path(out_dir, "prior_homo.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pri$hetersites,
            file = file.path(out_dir, "prior_heter.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pri$lambda1,
            file = file.path(out_dir, "prior_lambda1.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("done.\n")
