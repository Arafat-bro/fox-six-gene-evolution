#!/bin/bash
# Fo47 (GCA_013085055.1) is a required strain but was excluded from the
# automatic host-filtered shortlist because, as a nonpathogenic isolate,
# it has no associated host in NCBI's BioSample record. Add it manually
# to a separate "required additions" file.
set -euo pipefail

TABLE=results/selection/fo_genomes_table.tsv
ADDITIONS=results/selection/fo_genomes_manual_additions.tsv

head -1 "$TABLE" > "$ADDITIONS"
grep "GCA_013085055.1" "$TABLE" >> "$ADDITIONS"

echo "Manual additions file created:"
cat "$ADDITIONS" | column -t -s $'\t'
