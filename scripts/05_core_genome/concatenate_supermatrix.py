#!/usr/bin/env python3
"""
Part 9 (concatenation step).

Concatenates all successfully trimmed per-marker alignments into a
single core-genome supermatrix for IQ-TREE and creates a partition table.

Input:
    results/core_genome_markers/trimmed/*_trimmed.fasta

Outputs:
    results/core_genome_markers/core_genome_supermatrix.fasta
    results/core_genome_markers/core_genome_supermatrix.partitions.tsv
    results/core_genome_markers/concatenation_report.tsv

Usage:
    python3 scripts/05_core_genome/concatenate_supermatrix.py
"""

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


PROJECT_ROOT = Path.home() / "Foxy_SIX_Evolution"

TRIMMED_DIR = PROJECT_ROOT / "results/core_genome_markers/trimmed"
OUT_FASTA = PROJECT_ROOT / "results/core_genome_markers/core_genome_supermatrix.fasta"
OUT_PARTITION = PROJECT_ROOT / "results/core_genome_markers/core_genome_supermatrix.partitions.tsv"
OUT_REPORT = PROJECT_ROOT / "results/core_genome_markers/concatenation_report.tsv"


# Fixed genome order.
# This order is used consistently in the supermatrix and every downstream step.
GENOMES = [
    "GCA_003315725.1",
    "GCA_013085055.1",
    "GCA_003615085.1",
    "GCA_049306905.1",
    "GCA_053754965.2",
    "GCA_055944325.1",
]


def main():

    # ------------------------------------------------------------
    # Check input directory
    # ------------------------------------------------------------

    if not TRIMMED_DIR.exists():
        sys.exit(
            f"Trimmed alignment directory does not exist: {TRIMMED_DIR}\n"
            "Run align_trim.py first."
        )

    marker_files = sorted(TRIMMED_DIR.glob("*_trimmed.fasta"))

    if not marker_files:
        sys.exit(
            f"No trimmed alignments found in {TRIMMED_DIR}\n"
            "Run align_trim.py first."
        )

    print(f"Found {len(marker_files)} trimmed marker alignments.")
    print(f"Expected genomes per marker: {len(GENOMES)}")


    # ------------------------------------------------------------
    # Prepare output directory
    # ------------------------------------------------------------

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)


    # ------------------------------------------------------------
    # Storage
    # ------------------------------------------------------------

    concatenated = {acc: [] for acc in GENOMES}

    partition_rows = []
    report_rows = []

    current_pos = 1

    n_used = 0
    n_skipped = 0


    # ------------------------------------------------------------
    # Process every marker
    # ------------------------------------------------------------

    for marker_file in marker_files:

        marker_id = marker_file.stem.replace("_trimmed", "")

        records = {
            rec.id: str(rec.seq)
            for rec in SeqIO.parse(marker_file, "fasta")
        }

        record_ids = set(records)
        expected_ids = set(GENOMES)


        # --------------------------------------------------------
        # Check genome IDs
        # --------------------------------------------------------

        if record_ids != expected_ids:

            missing = sorted(expected_ids - record_ids)
            extra = sorted(record_ids - expected_ids)

            details = []

            if missing:
                details.append(
                    "missing genomes: " + ",".join(missing)
                )

            if extra:
                details.append(
                    "unexpected genomes: " + ",".join(extra)
                )

            report_rows.append(
                (marker_id, "SKIPPED", "; ".join(details))
            )

            n_skipped += 1
            continue


        # --------------------------------------------------------
        # Check alignment lengths
        # --------------------------------------------------------

        lengths = {
            len(records[acc])
            for acc in GENOMES
        }

        if len(lengths) != 1:

            report_rows.append(
                (
                    marker_id,
                    "SKIPPED",
                    f"unequal lengths: {sorted(lengths)}"
                )
            )

            n_skipped += 1
            continue


        marker_len = lengths.pop()


        # --------------------------------------------------------
        # Reject empty alignments
        # --------------------------------------------------------

        if marker_len == 0:

            report_rows.append(
                (marker_id, "SKIPPED", "zero-length alignment")
            )

            n_skipped += 1
            continue


        # --------------------------------------------------------
        # Add marker sequences to supermatrix
        # --------------------------------------------------------

        for acc in GENOMES:
            concatenated[acc].append(records[acc])


        # --------------------------------------------------------
        # Record partition coordinates
        # --------------------------------------------------------

        start = current_pos
        end = current_pos + marker_len - 1

        partition_rows.append(
            (marker_id, start, end)
        )

        report_rows.append(
            (
                marker_id,
                "USED",
                f"length={marker_len}"
            )
        )

        current_pos = end + 1
        n_used += 1


    # ------------------------------------------------------------
    # Check that at least one marker survived
    # ------------------------------------------------------------

    if n_used == 0:
        sys.exit(
            "No markers survived concatenation checks. "
            "Nothing to write."
        )


    # ------------------------------------------------------------
    # Build final supermatrix
    # ------------------------------------------------------------

    out_records = []

    for acc in GENOMES:

        full_seq = "".join(concatenated[acc])

        out_records.append(
            SeqRecord(
                Seq(full_seq),
                id=acc,
                description=""
            )
        )

    SeqIO.write(
        out_records,
        OUT_FASTA,
        "fasta"
    )


    # ------------------------------------------------------------
    # Write partition table
    # ------------------------------------------------------------

    with open(OUT_PARTITION, "w") as f:

        f.write("marker_id\tstart\tend\n")

        for marker_id, start, end in partition_rows:

            f.write(
                f"{marker_id}\t{start}\t{end}\n"
            )


    # ------------------------------------------------------------
    # Write concatenation report
    # ------------------------------------------------------------

    with open(OUT_REPORT, "w") as f:

        f.write(
            "marker_id\tstatus\tdetail\n"
        )

        for marker_id, status, detail in report_rows:

            f.write(
                f"{marker_id}\t{status}\t{detail}\n"
            )


    # ------------------------------------------------------------
    # Final checks
    # ------------------------------------------------------------

    final_length = current_pos - 1

    lengths_out = {
        acc: len("".join(seqs))
        for acc, seqs in concatenated.items()
    }


    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print()
    print("============================================================")
    print("CONCATENATION COMPLETE")
    print("============================================================")
    print(f"Markers found:          {len(marker_files)}")
    print(f"Markers used:           {n_used}")
    print(f"Markers skipped:        {n_skipped}")
    print(f"Genomes:                {len(GENOMES)}")
    print(f"Final supermatrix:      {final_length} columns")
    print()
    print(f"Supermatrix:")
    print(f"  {OUT_FASTA}")
    print()
    print(f"Partition table:")
    print(f"  {OUT_PARTITION}")
    print()
    print(f"Report:")
    print(f"  {OUT_REPORT}")
    print("============================================================")


    # ------------------------------------------------------------
    # Verify equal final sequence lengths
    # ------------------------------------------------------------

    if len(set(lengths_out.values())) != 1:

        print()
        print(
            "WARNING: genomes ended up with different total lengths."
        )

        for acc, length in lengths_out.items():
            print(f"  {acc}: {length}")

        sys.exit(
            "ERROR: final supermatrix is not column-consistent."
        )

    else:

        print()
        print(
            f"All {len(GENOMES)} genomes confirmed at equal length "
            f"({final_length} columns)."
        )

        print(
            "Concatenation is column-consistent."
        )


if __name__ == "__main__":
    main()
