#!/bin/bash
# Build the final genome metadata table required by the assignment:
# strain/isolate, forma specialis, host, geographic origin, accession,
# assembly level, genome size, contigs, N50, sequencing technology.
set -euo pipefail

LOOKUP=results/selection/forma_specialis_lookup.tsv
SELECTION=results/selection/fo_genomes_final_selection.tsv
OUT=results/genome_metadata.tsv

awk -F'\t' '
  NR==FNR { fs[$1]=$2; next }
  FNR==1 {
    print "Strain_Isolate\tForma_specialis\tHost\tGeographic_origin\tAssembly_Accession\tAssembly_Level\tGenome_Size_bp\tNum_Contigs\tContig_N50\tSequencing_Technology"
    next
  }
  {
    strain = ($3!="" ? $3 : ($4!="" ? $4 : ($6!="" ? $6 : $7)))
    accession = $1
    formaspecialis = (accession in fs ? fs[accession] : "")
    host = $8
    geo = $9
    if (host=="missing" || host=="not collected") host=""
    if (geo=="missing" || geo=="not collected") geo=""
    print strain"\t"formaspecialis"\t"host"\t"geo"\t"accession"\t"$11"\t"$13"\t"$14"\t"$15"\t"$12
  }
' "$LOOKUP" "$SELECTION" > "$OUT"

echo "Metadata table built. Row count (excluding header):"
tail -n +2 "$OUT" | wc -l

echo ""
column -t -s $'\t' "$OUT" | less -S
