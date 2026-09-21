#!/bin/bash
# Quick checks on the candidate genome pool:
# 1) How many genomes are Chromosome-level or better (our quality filter)?
# 2) Locate Fol4287 and Fo47 specifically (required strains).
set -euo pipefail

echo "=== Breakdown by Assembly Level ==="
cut -f4 results/selection/fo_genomes_table.tsv | sort | uniq -c

echo ""
echo "=== Searching for Fol4287 ==="
grep -i "4287" results/selection/fo_genomes_table.tsv || echo "Not found by name match — will need direct accession lookup."

echo ""
echo "=== Searching for Fo47 ==="
grep -iE "Fo47|Fo 47" results/selection/fo_genomes_table.tsv || echo "Not found by name match — will need direct accession lookup."
