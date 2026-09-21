#!/bin/bash
# Flatten downloaded genome FASTA files from NCBI's nested per-accession
# folder structure into a single flat directory, named by accession,
# for simpler use in downstream looped scripts.
set -euo pipefail

mkdir -p data/genomes/raw/flat

find data/genomes/raw/ncbi_dataset/data -name "*.fna" | while read -r filepath; do
  acc=$(basename "$(dirname "$filepath")")
  cp "$filepath" "data/genomes/raw/flat/${acc}.fasta"
done

echo "Flattened genome files:"
ls data/genomes/raw/flat | wc -l
