#!/bin/bash
# Trim the strawberry group to one representative per distinct growing
# region: retain Japan, Australia, and the highest-N50 California isolate;
# remove the two lower-N50 California isolates from the same region.
set -euo pipefail

FINAL=results/selection/fo_genomes_final_selection.tsv
TMP=results/selection/fo_genomes_final_selection_tmp.tsv

EXCLUDE_ACCESSIONS="GCA_016166095.2|GCA_016170085.2"

head -1 "$FINAL" > "$TMP"
tail -n +2 "$FINAL" | grep -vE "$EXCLUDE_ACCESSIONS" >> "$TMP"
mv "$TMP" "$FINAL"

echo "Updated final selection size:"
tail -n +2 "$FINAL" | wc -l
