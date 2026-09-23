"""
Apply manual scaffold-duplication exclusions that are NOT caught by either
automated dedup mechanism, on top of homolog_classification_final.tsv.

This is a distinct exclusion category from "excluded_pseudogene" (assigned
by 02_assign_final_status.py based on translation evidence) and from the
genuine-bucket-only manual dedup calls in
03_six_homolog_search/06_apply_manual_dedup_calls.py. Loci excluded here are
scaffold/haplotig duplicates discovered in the BORDERLINE bucket, which
neither automated dedup script ever inspects:
  - 05_dedupe_scaffold_duplicates.py only compares loci within the genuine
    bucket.
  - 06_apply_manual_dedup_calls.py reads only from
    genuine_homologs_deduped.tsv, so a manual exclusion added there for a
    borderline-bucket locus would silently have no effect.
This script is the correct place to apply such exclusions, since it acts on
the final classification table covering both genuine and borderline loci
together.

SIX1 / GCA_016164145.2: copy2 (WILW02000022.1:71194-72018, plus strand)
    excluded as a haplotig duplicate of copy3 (CM040507.1:2449316-2450140,
    minus strand). Confirmed by direct extraction and orientation-corrected
    comparison: sequence bodies byte-identical across the full ~825 bp
    locus (only FASTA headers differed). Both copies sit in the borderline
    identity bucket (64.06% identity, identical BLAST stats), which is
    exactly why neither automated dedup script flagged the pair. See
    docs/decisions/manual_duplication_calls.md for the full writeup.

Input:  results/six_homologs/classified/homolog_classification_final.tsv
Output: results/six_homologs/classified/homolog_classification_final.tsv
        (patched in place)
"""
import csv
from pathlib import Path

CLASS_FILE = Path("results/six_homologs/classified/homolog_classification_final.tsv")

MANUALLY_EXCLUDED_DUPLICATE_COPY_IDS = {
    "SIX1_GCA_016164145.2_copy2",  # haplotig duplicate of copy3 (CM040507.1);
                                     # missed by both automated dedup scripts
                                     # because both copies fell in the
                                     # borderline bucket. See
                                     # docs/decisions/manual_duplication_calls.md
}

with open(CLASS_FILE) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
    fieldnames = list(rows[0].keys())

n_changed = 0
for row in rows:
    if row["copy_id"] in MANUALLY_EXCLUDED_DUPLICATE_COPY_IDS and row["final_status"] == "retained":
        row["final_status"] = "excluded_manual_duplicate"
        n_changed += 1

with open(CLASS_FILE, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"Loci re-flagged as excluded_manual_duplicate: {n_changed}")
for cid in MANUALLY_EXCLUDED_DUPLICATE_COPY_IDS:
    match = [r for r in rows if r["copy_id"] == cid]
    if match:
        print(f"  {cid}: final_status = {match[0]['final_status']}")
    else:
        print(f"  {cid}: NOT FOUND in classification table -- check copy_id spelling")
