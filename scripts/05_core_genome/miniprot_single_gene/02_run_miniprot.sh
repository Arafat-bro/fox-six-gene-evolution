#!/bin/bash
mkdir -p results/core/gff
for g in data/genomes/raw/flat/*.fasta; do
  n=$(basename $g .fasta)
  [ -s results/core/gff/$n.gff ] && continue
  echo "start $n $(date)"
  miniprot -t 8 --gff -N 0 $g results/core/queries.faa > results/core/gff/$n.gff.tmp \
    && mv results/core/gff/$n.gff.tmp results/core/gff/$n.gff
  echo "done $n $(date)"
done