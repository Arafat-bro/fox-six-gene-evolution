#!/bin/bash
# Sequential BUSCO runner — one genome at a time, full 8 cores each.
# Rationale: after two full WSL2 restarts wiped 4-genome parallel runs,
# running sequentially means a crash only costs the ONE in-progress genome,
# not all of them, and each genome checkpoints to git as soon as it finishes.

set -uo pipefail  # NOT -e: we want to continue to the next genome even if one fails

GENOMES=(
  "GCA_003315725.1"
  "GCA_013085055.1"
  "GCA_003615085.1"
  "GCA_053754965.2"
  "GCA_055944325.1"
  "GCA_049306905.1"
)

LINEAGE_PATH="$HOME/busco_downloads/lineages/fungi_odb12.2"
GENOME_DIR="data/genomes/raw/flat"   # adjust if your flattened genome FASTAs live elsewhere
OUT_DIR="results/busco"

mkdir -p "$OUT_DIR" logs

for acc in "${GENOMES[@]}"; do
  echo "=========================================="
  echo "Starting $acc at $(date)"
  echo "=========================================="

  # Skip if already complete (idempotent — safe to re-run after a crash)
  if [ -f "$OUT_DIR/$acc/run_fungi_odb12.2/full_table.tsv" ]; then
    echo "$acc already has full_table.tsv — skipping."
    continue
  fi

  busco -i "$GENOME_DIR/${acc}.fasta" \
    -o "$acc" \
    --out_path "$OUT_DIR" \
    -l "$LINEAGE_PATH" \
    -m genome \
    -c 8 \
    --offline \
    -f \
    > "logs/${acc}.log" 2>&1

  status=$?

  if [ $status -eq 0 ] && [ -f "$OUT_DIR/$acc/run_fungi_odb12.2/full_table.tsv" ]; then
    echo "$acc finished successfully at $(date)"
    # Checkpoint: commit just the small summary output, not raw sequences
    git add "$OUT_DIR/$acc/run_fungi_odb12.2/full_table.tsv" \
            "$OUT_DIR/$acc/run_fungi_odb12.2/short_summary"*.txt \
            "logs/${acc}.log" 2>/dev/null
    git commit -m "BUSCO complete: $acc (fungi_odb12.2)" 2>/dev/null
    git push origin main 2>/dev/null
  else
    echo "$acc FAILED or produced no full_table.tsv — check logs/${acc}.log" >&2
    echo "Continuing to next genome anyway."
  fi
done

echo "=========================================="
echo "All genomes processed at $(date)"
echo "=========================================="
