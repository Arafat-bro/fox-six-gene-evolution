#!/usr/bin/env python3
"""
Classify each merged candidate SIX-gene locus as a genuine homolog,
a borderline case requiring manual review, or a rejected spurious hit,
using gene-specific identity and coverage thresholds.

See docs/decisions/homolog_identity_thresholds.md for the thresholds
used and their justification.

This script does NOT classify complete/truncated/pseudogenized/
assembly-boundary status -- that requires extracting and translating the
actual underlying sequence at each locus, which is a separate step.

Input:  results/six_homologs/merged/all_genes_merged_loci.tsv
Output: results/six_homologs/filtered/genuine_homologs.tsv
        results/six_homologs/filtered/borderline_for_manual_review.tsv
        results/six_homologs/filtered/rejected_below_threshold.tsv
"""

import csv
from pathlib import Path

IN_FILE = Path("results/six_homologs/merged/all_genes_merged_loci.tsv")
OUT_DIR = Path("results/six_homologs/filtered")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FIVE_GENES = {"SIX3", "SIX5", "SIX7", "SIX10", "SIX12"}

THRESHOLDS = {
    # gene_group: (identity_floor, identity_confident, coverage_floor_or_None)
    "five_genes": (70.0, 80.0, 70.0),
    "SIX1": (55.0, 65.0, None),
}


def thresholds_for(gene):
    return THRESHOLDS["SIX1"] if gene == "SIX1" else THRESHOLDS["five_genes"]


def classify(row):
    gene = row["gene"]
    identity = float(row["weighted_pident"])
    coverage = float(row["query_coverage_pct"])
    floor, confident, cov_floor = thresholds_for(gene)

    if identity < floor:
        return "rejected_below_identity_floor", "identity below gene-specific floor"

    if identity < confident:
        return "borderline_identity_manual_review", "identity between floor and confident threshold"

    if cov_floor is not None and coverage < cov_floor:
        return "genuine_low_coverage_candidate_truncation", \
               "confident identity but coverage below floor; possible truncation/assembly boundary, needs sequence check"

    return "genuine_homolog", ""


with open(IN_FILE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)

fieldnames = list(rows[0].keys()) + ["status", "status_note"] if rows else []

buckets = {
    "genuine_homolog": [],
    "genuine_low_coverage_candidate_truncation": [],
    "borderline_identity_manual_review": [],
    "rejected_below_identity_floor": [],
}

for row in rows:
    status, note = classify(row)
    row["status"] = status
    row["status_note"] = note
    buckets[status].append(row)

genuine_rows = buckets["genuine_homolog"] + buckets["genuine_low_coverage_candidate_truncation"]

with open(OUT_DIR / "genuine_homologs.tsv", "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(genuine_rows)

with open(OUT_DIR / "borderline_for_manual_review.tsv", "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(buckets["borderline_identity_manual_review"])

with open(OUT_DIR / "rejected_below_threshold.tsv", "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(buckets["rejected_below_identity_floor"])

print(f"Total candidate loci processed: {len(rows)}")
for status, bucket_rows in buckets.items():
    print(f"  {status}: {len(bucket_rows)}")

print("\nGenuine homologs (incl. low-coverage candidates) by gene:")
from collections import Counter
gene_counts = Counter(r["gene"] for r in genuine_rows)
for gene in ["SIX1", "SIX3", "SIX5", "SIX7", "SIX10", "SIX12"]:
    print(f"  {gene}: {gene_counts.get(gene, 0)}")
