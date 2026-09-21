#!/usr/bin/env python3
"""
Merge tblastn HSPs into candidate SIX-gene homolog loci per genome.

For spliced reference genes (SIX5, SIX10, SIX12), a single genomic copy is
frequently returned by tblastn as multiple separate HSPs, one per exon,
because BLAST has no knowledge of intron/exon structure. This script groups
HSPs that likely originate from a single underlying gene copy, based on
genomic proximity and non-overlapping query coverage, and computes an
overall query coverage and percent identity for each resulting candidate
locus (copy).

Merge rule (see docs/decisions/hsp_merging_criteria.md):
Two HSPs on the same scaffold and strand are merged into one candidate copy
only if (a) the genomic gap between them is <= 3000 bp AND (b) their query
(protein) ranges overlap by no more than 15 amino acids. HSPs that are
genomically close but whose query ranges overlap beyond this tolerance are
NOT merged, and are instead written to a flagged file for manual review,
since this pattern is also consistent with two independent nearby gene
copies rather than exon fragments of one gene.

Input:  results/six_homologs/raw/<GENE>_vs_<accession>.tsv
        data/references/processed/<GENE>_protein.fasta   (for query length)
Output: results/six_homologs/merged/<GENE>_vs_<accession>.merged.tsv
        results/six_homologs/merged/all_genes_merged_loci.tsv
        results/six_homologs/merged/flagged_for_manual_review.tsv
"""

import csv
from pathlib import Path
from Bio import SeqIO

RAW_DIR = Path("results/six_homologs/raw")
OUT_DIR = Path("results/six_homologs/merged")
QUERY_DIR = Path("data/references/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)

MAX_MERGE_GAP_BP = 3000
QUERY_OVERLAP_TOLERANCE_AA = 15

COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "sstrand"]


def load_query_length(gene, cache={}):
    if gene not in cache:
        fasta = QUERY_DIR / f"{gene}_protein.fasta"
        record = next(SeqIO.parse(fasta, "fasta"))
        cache[gene] = len(record.seq)
    return cache[gene]


def read_hsps(tsv_path):
    hsps = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, fieldnames=COLUMNS, delimiter="\t")
        for row in reader:
            row["pident"] = float(row["pident"])
            row["length"] = int(row["length"])
            row["qstart"] = int(row["qstart"])
            row["qend"] = int(row["qend"])
            row["sstart"] = int(row["sstart"])
            row["send"] = int(row["send"])
            row["evalue"] = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            hsps.append(row)
    return hsps


def genomic_span(hsp):
    return min(hsp["sstart"], hsp["send"]), max(hsp["sstart"], hsp["send"])


def query_span(hsp):
    return min(hsp["qstart"], hsp["qend"]), max(hsp["qstart"], hsp["qend"])


def overlap_len(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return max(0, hi - lo + 1)


def cluster_hsps(hsps, gene, accession, flagged):
    groups = {}
    for h in hsps:
        groups.setdefault((h["sseqid"], h["sstrand"]), []).append(h)

    clusters = []
    for (sseqid, sstrand), group in groups.items():
        group_sorted = sorted(group, key=lambda h: genomic_span(h)[0])
        current = [group_sorted[0]]
        current_hi = genomic_span(group_sorted[0])[1]
        for h in group_sorted[1:]:
            lo, hi = genomic_span(h)
            gap = lo - current_hi
            max_overlap = max(overlap_len(query_span(h), query_span(c)) for c in current)
            if gap <= MAX_MERGE_GAP_BP and max_overlap <= QUERY_OVERLAP_TOLERANCE_AA:
                current.append(h)
                current_hi = max(current_hi, hi)
            else:
                clusters.append((sseqid, sstrand, current))
                if gap <= MAX_MERGE_GAP_BP and max_overlap > QUERY_OVERLAP_TOLERANCE_AA:
                    flagged.append({
                        "gene": gene, "accession": accession,
                        "sseqid": sseqid, "sstrand": sstrand,
                        "genomic_gap_bp": gap, "query_overlap_aa": max_overlap,
                        "reason": "genomically close but query ranges overlap beyond tolerance "
                                  "(possible independent nearby copy, not an exon fragment)",
                    })
                current = [h]
                current_hi = hi
        clusters.append((sseqid, sstrand, current))
    return clusters


def union_query_coverage(hsps, query_len):
    intervals = sorted(query_span(h) for h in hsps)
    merged = []
    for lo, hi in intervals:
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    covered = sum(hi - lo + 1 for lo, hi in merged)
    return round(covered / query_len * 100, 1)


def weighted_identity(hsps):
    total_len = sum(h["length"] for h in hsps)
    return round(sum(h["pident"] * h["length"] for h in hsps) / total_len, 2) if total_len else 0.0


HEADER = ["copy_id", "gene", "accession", "sseqid", "sstrand", "n_hsps_merged",
          "genomic_start", "genomic_end", "query_coverage_pct", "weighted_pident",
          "best_evalue", "sum_bitscore"]

combined_rows = []
flagged_rows = []

for raw_file in sorted(RAW_DIR.glob("*.tsv")):
    gene, _, accession = raw_file.stem.partition("_vs_")
    query_len = load_query_length(gene)
    hsps = read_hsps(raw_file)
    out_path = OUT_DIR / f"{raw_file.stem}.merged.tsv"

    if not hsps:
        with open(out_path, "w") as fh:
            fh.write("\t".join(HEADER) + "\n")
        continue

    clusters = cluster_hsps(hsps, gene, accession, flagged_rows)
    clusters.sort(key=lambda c: min(genomic_span(h)[0] for h in c[2]))

    rows = []
    for i, (sseqid, sstrand, group) in enumerate(clusters, start=1):
        lo = min(genomic_span(h)[0] for h in group)
        hi = max(genomic_span(h)[1] for h in group)
        rows.append({
            "copy_id": f"{gene}_{accession}_copy{i}", "gene": gene, "accession": accession,
            "sseqid": sseqid, "sstrand": sstrand, "n_hsps_merged": len(group),
            "genomic_start": lo, "genomic_end": hi,
            "query_coverage_pct": union_query_coverage(group, query_len),
            "weighted_pident": weighted_identity(group),
            "best_evalue": min(h["evalue"] for h in group),
            "sum_bitscore": round(sum(h["bitscore"] for h in group), 1),
        })

    with open(out_path, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in HEADER) + "\n")
    combined_rows.extend(rows)

with open(OUT_DIR / "all_genes_merged_loci.tsv", "w") as fh:
    fh.write("\t".join(HEADER) + "\n")
    for r in combined_rows:
        fh.write("\t".join(str(r[c]) for c in HEADER) + "\n")

flag_header = ["gene", "accession", "sseqid", "sstrand", "genomic_gap_bp", "query_overlap_aa", "reason"]
with open(OUT_DIR / "flagged_for_manual_review.tsv", "w") as fh:
    fh.write("\t".join(flag_header) + "\n")
    for r in flagged_rows:
        fh.write("\t".join(str(r[c]) for c in flag_header) + "\n")

n_files = len(list(RAW_DIR.glob("*.tsv")))
print(f"Processed {n_files} raw BLAST result files.")
print(f"Total candidate loci (copies) identified: {len(combined_rows)}")
print(f"Flagged for manual review: {len(flagged_rows)}")
