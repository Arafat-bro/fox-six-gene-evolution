#!/bin/bash
# Filter candidate genomes to those that are Chromosome-level or Complete
# Genome (highest contiguity tiers) AND have a usable Host value, then
# save this as the working shortlist for manual selection.
set -euo pipefail

TABLE=results/selection/fo_genomes_table.tsv
OUT=results/selection/fo_genomes_shortlist.tsv

awk -F'\t' '
  NR==1 { print; next }
  ($11=="Chromosome" || $11=="Complete Genome") && $8!="" && $8!="missing" { print }
' "$TABLE" > "$OUT"

echo "Shortlisted genomes (Chromosome/Complete Genome level + known host):"
tail -n +2 "$OUT" | wc -l
