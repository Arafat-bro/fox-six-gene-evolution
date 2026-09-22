"""
Assign a single authoritative final_status column to every locus in
homolog_classification.tsv, so downstream steps (presence/absence matrix,
alignments, trees) can filter on one field instead of re-deriving inclusion
logic from raw classification columns each time.

Rules:
  - source_table == "borderline": retained only if classification is
    "complete" or "truncated" with no internal stops (i.e. the translation
    evidence supports treating it as a genuine divergent allele despite the
    low identity score that originally flagged it). Pseudogenized borderline
    loci are excluded, consistent with the SIX1_GCA_053754955.2_copy1 manual
    call (see docs/decisions/manual_duplication_calls.md and the walkthrough
    of that case).
  - source_table == "genuine": retained if classification is "complete" or
    "truncated"; excluded if "pseudogenized". Complete/truncated genuine loci
    keep their status regardless of assembly_boundary_flag -- the boundary
    flag is carried forward as a caveat column, not used to exclude here,
    since a boundary-flagged locus may still be a real (if unconfirmed)
    detection rather than a false one.

Input:  results/six_homologs/classified/homolog_classification.tsv
Output: results/six_homologs/classified/homolog_classification_final.tsv
"""
import csv
from pathlib import Path

IN_FILE = Path("results/six_homologs/classified/homolog_classification.tsv")
OUT_FILE = Path("results/six_homologs/classified/homolog_classification_final.tsv")


def final_status(row):
    if row["classification"] == "pseudogenized":
        return "excluded_pseudogene"
    return "retained"


with open(IN_FILE) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
    fieldnames = list(rows[0].keys()) + ["final_status"]

for row in rows:
    row["final_status"] = final_status(row)

with open(OUT_FILE, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"Total loci: {len(rows)}")
print(f"  retained: {sum(1 for r in rows if r['final_status'] == 'retained')}")
print(f"  excluded_pseudogene: {sum(1 for r in rows if r['final_status'] == 'excluded_pseudogene')}")
print()
print("By source_table x final_status:")
from collections import Counter
print(Counter((r["source_table"], r["final_status"]) for r in rows))
