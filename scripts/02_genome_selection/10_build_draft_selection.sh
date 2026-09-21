#!/bin/bash
# Assemble the working draft genome selection: start from the plant-host
# shortlist, remove genomes excluded for redundancy (flax, banana), and
# add the manually-included required strain (Fo47).
set -euo pipefail

SHORTLIST=results/selection/fo_genomes_shortlist_plants.tsv
ADDITIONS=results/selection/fo_genomes_manual_additions.tsv
OUT=results/selection/fo_genomes_draft_selection.tsv

EXCLUDE_ACCESSIONS="GCA_013423225.1|GCA_013423235.1|GCA_013423255.1|GCA_013423265.1|GCA_011316005.3|GCA_014282265.3"

head -1 "$SHORTLIST" > "$OUT"
tail -n +2 "$SHORTLIST" | grep -vE "$EXCLUDE_ACCESSIONS" >> "$OUT"
tail -n +2 "$ADDITIONS" >> "$OUT"

echo "Draft selection size:"
tail -n +2 "$OUT" | wc -l

echo ""
echo "=== Draft selection, sorted by host ==="
awk -F'\t' 'NR>1 { print $8"\t"$1"\t"$3"\t"$11"\t"$15 }' "$OUT" | sort | column -t -s $'\t'
