#!/usr/bin/env bash
# Run-time benchmark: Rust port vs the R reference on the bundled PAR-CLIP
# example, then on a synthesized 100k-row dataset for a representative
# real-data load. All commands run with seed=123, bg.method=lrt,
# bg.target=both, lambda1=auto, and thread=1 unless explicitly varied.
#
# To make the comparison fair, both the R and Rust runs are restricted to
# `estimate_inference_with_empirical` only, with no compute_prior or
# prior-aware inference in the timed region.

set -euo pipefail
cd "$(dirname "$0")/.."

PARCLIP=validation/parclip.tsv
LARGE=validation/parclip_100k.tsv
RUST=./target/release/mirage
mkdir -p validation/bench_out

cargo build --release >/dev/null 2>&1

if [[ ! -f "$LARGE" ]]; then
  python3 - <<'PY' > "$LARGE"
with open("validation/parclip.tsv") as f:
    header = f.readline()
    rows = f.readlines()
import sys
sys.stdout.write(header)
for k in range(100):
    for r in rows:
        cols = r.rstrip("\n").split("\t")
        cols[0] = f"{cols[0]}_rep{k}"
        sys.stdout.write("\t".join(cols) + "\n")
PY
fi

bench_R_empirical_only () {
  local input="$1"; local label="$2"
  /usr/bin/env time -f "${label}: %e s (user %U, sys %S, max-rss %M KB)" \
    Rscript -e "
suppressPackageStartupMessages(library(MIRAGE))
ct <- as.data.frame(data.table::fread('$input'))
t0 <- proc.time()
res <- with(ct, estimate_inference_with_empirical(
  cbind(pos, motif, type),
  treated_fixed_count, control_fixed_count, treated_depth, control_depth,
  delta=0.001, depth.cutoff=10, lambda1='auto', top.sites=NULL,
  bg.method='lrt', bg.target='both', highly.methyl.cutoff=0.95,
  seed=123, thread=1))
cat(sprintf('R inner time: %.3f s\n', (proc.time() - t0)[3]))
invisible(res)
" >/dev/null
}

bench_rust_empirical_only () {
  local input="$1"; local thread="$2"; local label="$3"; local out="$4"
  /usr/bin/env time -f "${label}: %e s (user %U, sys %S, max-rss %M KB)" \
    "$RUST" estimate-empirical \
      --counts "$input" \
      --delta 0.001 --depth-cutoff 10 --lambda1 auto \
      --bg-method lrt --bg-target both --highly-methyl-cutoff 0.95 \
      --seed 123 --thread "$thread" \
      --out-homo "validation/bench_out/${out}_homo.tsv" \
      --out-heter "validation/bench_out/${out}_heter.tsv" \
      --out-summary "validation/bench_out/${out}_summary.tsv" >/dev/null
}

echo "===== PAR-CLIP example (999 sites) ====="
bench_R_empirical_only "$PARCLIP" "R   t=1"
for t in 1 4 8; do
  bench_rust_empirical_only "$PARCLIP" "$t" "Rust t=$t" "small_t${t}"
done

echo
echo "===== Synthetic 100k sites (tiled PAR-CLIP) ====="
echo "Lines: $(wc -l < "$LARGE")"
bench_R_empirical_only "$LARGE" "R   t=1"
for t in 1 4 8; do
  bench_rust_empirical_only "$LARGE" "$t" "Rust t=$t" "large_t${t}"
done
