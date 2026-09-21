#!/bin/bash
# Build a nucleotide BLAST database for each of the 27 selected genomes.
# Input:  data/genomes/raw/flat/<accession>.fasta
# Output: BLAST db files alongside each input FASTA (<accession>.fasta.n*)
# Log:    logs/build_blast_databases.log

set -euo pipefail

GENOME_DIR="data/genomes/raw/flat"
LOG_FILE="logs/build_blast_databases.log"

> "$LOG_FILE"

count=0
for genome in "$GENOME_DIR"/*.fasta; do
  acc=$(basename "$genome" .fasta)
  echo "Building BLAST db for $acc ..." | tee -a "$LOG_FILE"
  makeblastdb -in "$genome" -dbtype nucl -parse_seqids >> "$LOG_FILE" 2>&1
  count=$((count + 1))
done

echo "Done. Built $count BLAST databases." | tee -a "$LOG_FILE"
