#!/bin/bash
# Download the actual genome assembly FASTA files for the final
# 27-genome selection.
set -euo pipefail

mkdir -p data/genomes/raw

datasets download genome accession \
  --inputfile scripts/02_genome_selection/final_accessions.txt \
  --include genome

unzip -o ncbi_dataset.zip -d data/genomes/raw/

echo "Download complete. Genome files:"
find data/genomes/raw -name "*.fna" | wc -l
