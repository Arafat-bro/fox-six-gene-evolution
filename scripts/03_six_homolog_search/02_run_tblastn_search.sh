#!/bin/bash
# Run tblastn for each of the six SIX reference proteins against each of the
# 27 selected genome BLAST databases.
# Input:  data/references/processed/<GENE>_protein.fasta   (6 files)
#         data/genomes/raw/flat/<accession>.fasta          (27 BLAST dbs, built earlier)
# Output: results/six_homologs/raw/<GENE>_vs_<accession>.tsv   (162 files)
# Log:    logs/run_tblastn_search.log

set -euo pipefail

QUERY_DIR="data/references/processed"
DB_DIR="data/genomes/raw/flat"
OUT_DIR="results/six_homologs/raw"
LOG_FILE="logs/run_tblastn_search.log"

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs sstrand"

mkdir -p "$OUT_DIR"
> "$LOG_FILE"

count=0
for query in "$QUERY_DIR"/*_protein.fasta; do
  gene=$(basename "$query" _protein.fasta)
  for genome in "$DB_DIR"/*.fasta; do
    acc=$(basename "$genome" .fasta)
    out_file="$OUT_DIR/${gene}_vs_${acc}.tsv"

    echo "Searching $gene vs $acc ..." | tee -a "$LOG_FILE"
    tblastn -query "$query" -db "$genome" -outfmt "$OUTFMT" -evalue 1e-10 > "$out_file" 2>> "$LOG_FILE"

    count=$((count + 1))
  done
done

echo "Done. Ran $count tblastn searches." | tee -a "$LOG_FILE"
