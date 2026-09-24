#!/usr/bin/env python3
"""
Filter and intersect BUSCO single-copy orthologs across all core-genome-tree genomes.

WHY THIS SCRIPT EXISTS (for methods section / defending to professor):
A marker gene is only usable for the concatenated supermatrix if it is:
  1. Complete AND single-copy in EVERY genome (not just present somewhere) —
     a gene missing or duplicated in even one genome can't go into a single
     concatenated alignment column-for-column across all genomes.
  2. Free of internal stop codons in EVERY genome — an internal stop is a
     miniprot gene-prediction artifact (BUSCO flags these directly), not a
     real biological truncation. Including it would inject broken/garbage
     sequence into the alignment and downstream tree.

This is a straightforward set-intersection problem, done in code rather than
by hand so it's reproducible and auditable: reviewers can rerun this script
against the same full_table.tsv files and get the identical marker list.
"""

import csv
import sys
from pathlib import Path

# ---- CONFIG: your 6 core-genome-tree genomes ----
GENOMES = [
    "GCA_003315725.1",   # Fol4287
    "GCA_013085055.1",   # Fo47
    "GCA_003615085.1",   # FoC_Fus2
    "GCA_049306905.1",
    "GCA_053754965.2",
    "GCA_055944325.1",
]
LINEAGE_RUN = "run_hypocreales_odb12.2"
BUSCO_BASE = Path("results/busco")
OUT_DIR = Path("results/core_genome_markers")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_full_table(path):
    """Return set of BUSCO IDs marked 'Complete' (i.e. single-copy complete;
    BUSCO labels duplicated hits as 'Duplicated', a different status, so this
    filter alone already excludes duplicates)."""
    complete_ids = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            busco_id, status = fields[0], fields[1]
            if status == "Complete":
                complete_ids.add(busco_id)
    return complete_ids


def find_stop_codon_ids(genome_dir, complete_ids):
    """Scan single-copy sequence FASTAs for internal stop codons ('*' before
    the final residue). Returns the subset of complete_ids that are
    contaminated and must be excluded."""
    seq_dir = genome_dir / LINEAGE_RUN / "busco_sequences" / "single_copy_busco_sequences"
    flagged = set()
    for busco_id in complete_ids:
        faa_path = seq_dir / f"{busco_id}.faa"
        if not faa_path.exists():
            # Shouldn't happen for a Complete single-copy hit, but don't
            # silently assume clean -- flag for manual look if missing.
            print(f"  WARNING: expected sequence file missing: {faa_path}", file=sys.stderr)
            flagged.add(busco_id)
            continue
        seq_lines = []
        with open(faa_path) as f:
            for line in f:
                if not line.startswith(">"):
                    seq_lines.append(line.strip())
        seq = "".join(seq_lines)
        # Internal stop = '*' anywhere except possibly the last character.
        if "*" in seq[:-1] or (seq.endswith("*") and seq.count("*") > 1):
            flagged.add(busco_id)
    return flagged


def main():
    per_genome_clean = {}
    summary_rows = []

    for acc in GENOMES:
        genome_dir = BUSCO_BASE / acc
        full_table = genome_dir / LINEAGE_RUN / "full_table.tsv"
        if not full_table.exists():
            print(f"FATAL: missing full_table.tsv for {acc} at {full_table}", file=sys.stderr)
            sys.exit(1)

        complete_ids = parse_full_table(full_table)
        stop_flagged = find_stop_codon_ids(genome_dir, complete_ids)
        clean_ids = complete_ids - stop_flagged
        per_genome_clean[acc] = clean_ids

        summary_rows.append({
            "genome": acc,
            "complete_single_copy": len(complete_ids),
            "internal_stop_flagged": len(stop_flagged),
            "clean_single_copy": len(clean_ids),
        })
        print(f"{acc}: {len(complete_ids)} complete/single-copy, "
              f"{len(stop_flagged)} flagged (internal stop), "
              f"{len(clean_ids)} clean")

    # Intersection across ALL genomes: a marker only qualifies for the
    # supermatrix if it's clean in every single genome.
    intersection = set.intersection(*per_genome_clean.values())
    print(f"\nFinal intersection across all {len(GENOMES)} genomes: "
          f"{len(intersection)} markers")

    # Write outputs
    marker_list_path = OUT_DIR / "core_genome_marker_ids.txt"
    with open(marker_list_path, "w") as f:
        for busco_id in sorted(intersection):
            f.write(busco_id + "\n")

    summary_path = OUT_DIR / "per_genome_filtering_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "genome", "complete_single_copy", "internal_stop_flagged", "clean_single_copy"
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nWrote marker ID list -> {marker_list_path}")
    print(f"Wrote per-genome summary -> {summary_path}")


if __name__ == "__main__":
    main()
