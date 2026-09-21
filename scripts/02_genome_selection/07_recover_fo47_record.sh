#!/bin/bash
# Fo47 is the required nonpathogenic control strain. It was excluded by
# the host-data filter because nonpathogenic isolates typically have no
# associated host in NCBI's BioSample record. Recover its full record
# from the unfiltered table for manual inclusion.
set -euo pipefail

TABLE=results/selection/fo_genomes_table.tsv

echo "=== All Fo47 records (any assembly level) ==="
awk -F'\t' 'NR==1 || /Fo47/' "$TABLE" | column -t -s $'\t'
