#!/bin/bash
# Sanity check: confirm every N50 value in the final metadata table
# matches the original NCBI table, using a single-pass awk join
# rather than a per-row shell loop (which is error-prone for this).
set -euo pipefail

ORIGINAL=results/selection/fo_genomes_table.tsv
FINAL=results/genome_metadata.tsv

awk -F'\t' '
  NR==FNR { n50[$1]=$15; next }
  FNR==1 { next }
  {
    acc=$5
    if (!(acc in n50)) { print "NOT FOUND in original table: "acc; next }
    if ($9 != n50[acc]) { print "MISMATCH for "acc": final="$9", original="n50[acc] }
  }
' "$ORIGINAL" "$FINAL"

echo "Check complete. No output above (besides this line) means all N50 values match."
