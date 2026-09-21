#!/bin/bash
# Finalize the genome selection: remove the redundant duplicate Brassica
# assembly (same isolate, lower-quality version) and the ambiguous
# banana "race 4" genome whose specific race could not be confirmed.
set -euo pipefail

DRAFT=results/selection/fo_genomes_draft_selection.tsv
FINAL=results/selection/fo_genomes_final_selection.tsv

EXCLUDE_ACCESSIONS="GCA_014154955.1|GCA_027920445.1"

head -1 "$DRAFT" > "$FINAL"
tail -n +2 "$DRAFT" | grep -vE "$EXCLUDE_ACCESSIONS" >> "$FINAL"

echo "Final selection size:"
tail -n +2 "$FINAL" | wc -l
