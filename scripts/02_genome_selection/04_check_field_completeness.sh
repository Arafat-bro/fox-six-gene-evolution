#!/bin/bash
# Check how many of the 885 candidate genomes have usable (non-blank,
# non-"missing") values for the metadata fields the assignment requires,
# to determine how much manual BioSample review will be needed.
set -euo pipefail

TABLE=results/selection/fo_genomes_table.tsv
TOTAL=$(($(wc -l < $TABLE) - 1))  # subtract 1 for header row

echo "Total candidate genomes: $TOTAL"
echo ""

# Column numbers, per the header: 8=Host, 9=Geographic location, 12=Sequencing Tech
for col_num in 8 9 12; do
  col_name=$(head -1 $TABLE | cut -f$col_num)
  count=$(tail -n +2 $TABLE | cut -f$col_num | grep -vic '^$\|^missing$')
  echo "$col_name: $count / $TOTAL have a usable value"
done
