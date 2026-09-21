#!/bin/bash
# Remove non-plant hosts from the shortlist. This project focuses on
# plant-pathogenic formae speciales; isolates from non-plant hosts
# (e.g., clinical/human isolates) fall outside its scope.
set -euo pipefail

SHORTLIST=results/selection/fo_genomes_shortlist.tsv
OUT=results/selection/fo_genomes_shortlist_plants.tsv

awk -F'\t' '
  NR==1 { print; next }
  $8 != "Homo sapiens" { print }
' "$SHORTLIST" > "$OUT"

echo "Removed non-plant host records. Remaining:"
tail -n +2 "$OUT" | wc -l
