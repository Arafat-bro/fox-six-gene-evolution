#!/bin/bash
# Prints FASTA header and CDS feature block for each reference gene,
# for manual inspection of completeness/intron structure.
# Findings are recorded in docs/decisions/reference_gene_verification.md
cd ../../data/references/raw

for acc in MK906592.1 MK906598.1 MK906607.1 GQ268954.1 MK906667.1 MW160867.1; do
  echo "=== $acc ==="
  head -1 ${acc}.fasta
  grep -A 3 "CDS" ${acc}.gb
  echo
done
