#!/bin/bash
# Downloads FASTA + GenBank flatfile for the six SIX gene reference accessions.
set -e
cd ../../data/references/raw

for acc in MK906592.1 MK906598.1 MK906607.1 GQ268954.1 MK906667.1 MW160867.1; do
  efetch -db nuccore -id $acc -format fasta > ${acc}.fasta
  efetch -db nuccore -id $acc -format gb > ${acc}.gb
done
