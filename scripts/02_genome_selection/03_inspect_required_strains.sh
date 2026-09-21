#!/bin/bash
# Display candidate records for the two required strains (Fol4287, Fo47)
# in clean, aligned columns for manual comparison and selection.
set -euo pipefail

echo "=== Fol4287 candidates ==="
awk -F'\t' 'NR==1 || /4287/' results/selection/fo_genomes_table.tsv | column -t -s $'\t'

echo ""
echo "=== Fo47 candidates ==="
awk -F'\t' 'NR==1 || /Fo47/' results/selection/fo_genomes_table.tsv | column -t -s $'\t'
