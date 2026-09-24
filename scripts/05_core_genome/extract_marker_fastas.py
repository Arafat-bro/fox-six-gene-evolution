#!/usr/bin/env python3
"""
Part 9 (extraction step) - Foxy_SIX_Evolution core-genome pipeline.

For every core-genome marker ID (the 3620-marker intersection), pull that
marker's single-copy BUSCO amino-acid sequence out of each of the 6 genomes'
BUSCO run directories, relabel the FASTA header to the genome's accession
(not BUSCO's internal marker-specific ID), and write one per-marker
multi-FASTA containing all 6 genomes' copies.

Usage:
    python3 extract_marker_fastas.py

Reads:
    results/core_genome_markers/core_genome_marker_ids.txt
    results/busco/<ACCESSION>/run_hypocreales_odb12.2/busco_sequences/
        single_copy_busco_sequences/<marker_id>.faa

Writes:
    results/core_genome_markers/per_marker_fasta/<marker_id>.fasta
    results/core_genome_markers/extraction_report.tsv
"""

import sys
from pathlib import Path
from Bio import SeqIO

PROJECT_ROOT = Path.home() / "Foxy_SIX_Evolution"
MARKER_ID_FILE = PROJECT_ROOT / "results/core_genome_markers/core_genome_marker_ids.txt"
BUSCO_RESULTS_DIR = PROJECT_ROOT / "results/busco"
OUT_DIR = PROJECT_ROOT / "results/core_genome_markers/per_marker_fasta"
REPORT_FILE = PROJECT_ROOT / "results/core_genome_markers/extraction_report.tsv"

# Explicit list, not a glob over results/busco/* - this guarantees a stray
# extra genome folder can never silently sneak into the tree. A 7th genome
# appearing partway through would be a very hard bug to notice later.
GENOMES = [
    "GCA_003315725.1",  # Fol4287, reference
    "GCA_013085055.1",  # Fo47, nonpathogenic control
    "GCA_003615085.1",  # FoC_Fus2, f.sp. cepae
    "GCA_049306905.1",
    "GCA_053754965.2",
    "GCA_055944325.1",
]

BUSCO_LINEAGE_DIR = "run_hypocreales_odb12.2"


def busco_faa_dir(accession: str) -> Path:
    return (
        BUSCO_RESULTS_DIR
        / accession
        / BUSCO_LINEAGE_DIR
        / "busco_sequences"
        / "single_copy_busco_sequences"
    )


def main():
    if not MARKER_ID_FILE.exists():
        sys.exit(f"Marker ID file not found: {MARKER_ID_FILE}")

    marker_ids = [
        line.strip() for line in MARKER_ID_FILE.read_text().splitlines() if line.strip()
    ]
    print(f"Loaded {len(marker_ids)} marker IDs from {MARKER_ID_FILE}")

    # Fail fast if a genome's BUSCO directory is missing entirely - better to
    # find out now than 3600 markers into a partially-broken run.
    for acc in GENOMES:
        d = busco_faa_dir(acc)
        if not d.is_dir():
            sys.exit(f"Missing single_copy_busco_sequences dir for {acc}: {d}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    n_complete = 0
    n_incomplete = 0
    incomplete_markers = []

    with open(REPORT_FILE, "w") as report:
        report.write("marker_id\tgenomes_found\tstatus\n")

        for marker_id in marker_ids:
            records = []
            missing_from = []

            for acc in GENOMES:
                faa_path = busco_faa_dir(acc) / f"{marker_id}.faa"
                if not faa_path.exists():
                    missing_from.append(acc)
                    continue
                seqs = list(SeqIO.parse(faa_path, "fasta"))
                if len(seqs) != 1:
                    # A marker file with 0 or >1 sequences means BUSCO's
                    # single-copy claim doesn't hold here after all -
                    # shouldn't happen post-intersection, but check rather
                    # than assume.
                    missing_from.append(f"{acc}(n={len(seqs)})")
                    continue
                rec = seqs[0]
                rec.id = acc
                rec.description = ""
                records.append(rec)

            if len(records) == len(GENOMES):
                out_path = OUT_DIR / f"{marker_id}.fasta"
                SeqIO.write(records, out_path, "fasta")
                report.write(f"{marker_id}\t{len(records)}\tOK\n")
                n_complete += 1
            else:
                report.write(
                    f"{marker_id}\t{len(records)}\tINCOMPLETE (missing: {','.join(missing_from)})\n"
                )
                n_incomplete += 1
                incomplete_markers.append(marker_id)

    print(f"Complete (6/6 genomes): {n_complete}")
    print(f"Incomplete (dropped): {n_incomplete}")
    if incomplete_markers:
        print("First few dropped markers:", incomplete_markers[:5])
    print(f"Full report: {REPORT_FILE}")

    if n_incomplete > 0:
        print(
            "\nNOTE: the intersection step was supposed to guarantee 6/6 "
            "presence for every marker. Any INCOMPLETE rows here mean "
            "something is inconsistent between core_genome_marker_ids.txt "
            "and the actual .faa files on disk - worth checking before "
            "proceeding, not just accepting a smaller final number."
        )


if __name__ == "__main__":
    main()
