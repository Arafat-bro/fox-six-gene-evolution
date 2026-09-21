#!/bin/bash
# Fetch all F. oxysporum genome assembly records from NCBI Datasets,
# including BioSample-derived metadata (host, geographic origin,
# sequencing technology), and convert to a flat TSV table.
set -euo pipefail

mkdir -p results/selection

datasets summary genome taxon "Fusarium oxysporum" \
  --assembly-source GenBank \
  --as-json-lines > results/selection/fo_genomes_raw.jsonl

dataformat tsv genome \
  --inputfile results/selection/fo_genomes_raw.jsonl \
  --fields accession,organism-name,organism-infraspecific-strain,organism-infraspecific-isolate,assminfo-biosample-accession,assminfo-biosample-strain,assminfo-biosample-isolate,assminfo-biosample-host,assminfo-biosample-geo-loc-name,assminfo-biosample-collection-date,assminfo-level,assminfo-sequencing-tech,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50 \
  > results/selection/fo_genomes_table.tsv

echo "Total genome records found:"
wc -l < results/selection/fo_genomes_table.tsv
