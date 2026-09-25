#!/usr/bin/env python3
"""
Part 9 (alignment + trimming step).

Runs MAFFT (--auto) then trimAl (-automated1) on every per-marker
multi-FASTA from extract_marker_fastas.py, parallelized across markers.

Why parallelize this step but not BUSCO: BUSCO's memory problem was
miniprot_align scanning a whole genome (~5.5GB/process). Here each job
aligns just 6 short protein sequences - trivially light, so running
several in parallel across 8 cores is safe.

Usage:
    python3 align_trim.py [--workers N]
"""

import argparse
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

PROJECT_ROOT = Path.home() / "Foxy_SIX_Evolution"
IN_DIR = PROJECT_ROOT / "results/core_genome_markers/per_marker_fasta"
ALIGNED_DIR = PROJECT_ROOT / "results/core_genome_markers/aligned"
TRIMMED_DIR = PROJECT_ROOT / "results/core_genome_markers/trimmed"
LOG_FILE = PROJECT_ROOT / "results/core_genome_markers/align_trim_failures.tsv"


def align_and_trim(marker_id: str):
    """Runs in a worker process. Returns (marker_id, status, detail)."""
    in_fasta = IN_DIR / f"{marker_id}.fasta"
    aligned_fasta = ALIGNED_DIR / f"{marker_id}_aligned.fasta"
    trimmed_fasta = TRIMMED_DIR / f"{marker_id}_trimmed.fasta"

    try:
        with open(aligned_fasta, "w") as out:
            result = subprocess.run(
                ["mafft", "--auto", "--quiet", str(in_fasta)],
                stdout=out, stderr=subprocess.PIPE, timeout=300,
            )
        if result.returncode != 0:
            return marker_id, "MAFFT_FAILED", result.stderr.decode()[:300]
    except subprocess.TimeoutExpired:
        return marker_id, "MAFFT_TIMEOUT", ""

    try:
        result = subprocess.run(
            ["trimal", "-in", str(aligned_fasta), "-out", str(trimmed_fasta), "-automated1"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120,
        )
        if result.returncode != 0:
            return marker_id, "TRIMAL_FAILED", result.stderr.decode()[:300]
    except subprocess.TimeoutExpired:
        return marker_id, "TRIMAL_TIMEOUT", ""

    # trimAl can, in principle, trim a short/gappy/variable locus down to
    # zero columns - check for this rather than silently carrying an empty
    # alignment into the supermatrix.
    if trimmed_fasta.stat().st_size == 0:
        return marker_id, "TRIMMED_TO_EMPTY", ""

    return marker_id, "OK", ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6,
                         help="Parallel workers (default 6 of 8 cores - leaves headroom for WSL/OS).")
    args = parser.parse_args()

    ALIGNED_DIR.mkdir(parents=True, exist_ok=True)
    TRIMMED_DIR.mkdir(parents=True, exist_ok=True)

    marker_ids = sorted(p.stem for p in IN_DIR.glob("*.fasta"))
    if not marker_ids:
        sys.exit(f"No per-marker FASTAs found in {IN_DIR} - run extract_marker_fastas.py first.")
    print(f"Aligning + trimming {len(marker_ids)} markers with {args.workers} workers...")

    n_ok = 0
    failures = []

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(align_and_trim, m): m for m in marker_ids}
        for i, fut in enumerate(as_completed(futures), 1):
            marker_id, status, detail = fut.result()
            if status == "OK":
                n_ok += 1
            else:
                failures.append((marker_id, status, detail))
            if i % 200 == 0 or i == len(marker_ids):
                print(f"  {i}/{len(marker_ids)} done ({n_ok} OK, {len(failures)} failed)")

    with open(LOG_FILE, "w") as f:
        f.write("marker_id\tstatus\tdetail\n")
        for marker_id, status, detail in failures:
            f.write(f"{marker_id}\t{status}\t{detail}\n")

    print(f"\nDone. OK: {n_ok}  Failed/dropped: {len(failures)}")
    if failures:
        print(f"Failure breakdown logged to {LOG_FILE}")
        from collections import Counter
        counts = Counter(s for _, s, _ in failures)
        for status, count in counts.most_common():
            print(f"  {status}: {count}")
        print(
            "\nMarkers that failed or trimmed to empty are simply excluded "
            "from the supermatrix - the final marker count will land "
            "slightly below 3620, and that number + reason belongs in your "
            "methods section, not silently absorbed."
        )


if __name__ == "__main__":
    main()
