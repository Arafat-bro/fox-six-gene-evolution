#!/bin/bash
# Display the shortlist grouped by host, and confirm whether the two
# required strains (Fol4287, Fo47) are present in it.
set -euo pipefail

SHORTLIST=results/selection/fo_genomes_shortlist.tsv

echo "=== Shortlist sorted by Host (accession, strain, host, level, N50) ==="
awk -F'\t' 'NR>1 { print $8"\t"$1"\t"$3"\t"$11"\t"$15 }' "$SHORTLIST" | sort | column -t -s $'\t'

echo ""
echo "=== Host value counts (how many genomes per host) ==="
awk -F'\t' 'NR>1 { print $8 }' "$SHORTLIST" | sort | uniq -c | sort -rn

echo ""
echo "=== Checking for our chosen Fol4287 accession (GCA_003315725.1) ==="
grep "GCA_003315725.1" "$SHORTLIST" || echo "NOT in shortlist — will need manual override."

echo ""
echo "=== Checking for our chosen Fo47 accession (GCA_013085055.1) ==="
grep "GCA_013085055.1" "$SHORTLIST" || echo "NOT in shortlist — will need manual override."
