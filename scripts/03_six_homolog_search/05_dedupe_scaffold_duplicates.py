#!/usr/bin/env python3
"""
Check every (gene, genome) pair with more than one surviving locus for
assembly-duplicated scaffolds: cases where two loci represent the same
underlying genomic sequence assembled twice (identical or reverse-
complement-identical sequence), rather than a true second gene copy.

Input:  results/six_homologs/filtered/genuine_homologs.tsv
        data/genomes/raw/flat/<accession>.fasta
Output: results/six_homologs/filtered/genuine_homologs_deduped.tsv
        docs/decisions/scaffold_duplication_check_auto.tsv
            (log of every multi-locus pair checked and the verdict)
"""

import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

IN_FILE = Path("results/six_homologs/filtered/genuine_homologs.tsv")
GENOME_DIR = Path("data/genomes/raw/flat")
OUT_FILE = Path("results/six_homologs/filtered/genuine_homologs_deduped.tsv")
LOG_FILE = Path("docs/decisions/scaffold_duplication_check_auto.tsv")


def clean_seqid(sseqid):
    # sseqid looks like gb|WILW02000022.1| -- extract the accession part
    parts = sseqid.strip("|").split("|")
    return parts[1] if len(parts) > 1 else sseqid


def get_sequence(genome_records, sseqid, start, end):
    acc = clean_seqid(sseqid)
    if acc not in genome_records:
        return None
    lo, hi = min(start, end), max(start, end)
    return str(genome_records[acc].seq[lo - 1:hi])


def sequences_match(seq_a, seq_b):
    if seq_a is None or seq_b is None:
        return False
    if seq_a == seq_b:
        return "identical"
    if seq_a == str(Seq(seq_b).reverse_complement()):
        return "reverse_complement"
    return False


with open(IN_FILE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames

groups = {}
for row in rows:
    key = (row["gene"], row["accession"])
    groups.setdefault(key, []).append(row)

genome_cache = {}
log_rows = []
excluded_copy_ids = set()

for (gene, accession), group in groups.items():
    if len(group) < 2:
        continue

    if accession not in genome_cache:
        fasta_path = GENOME_DIR / f"{accession}.fasta"
        genome_cache[accession] = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    genome_records = genome_cache[accession]

    sequences = {}
    for row in group:
        seq = get_sequence(genome_records, row["sseqid"], int(row["genomic_start"]), int(row["genomic_end"]))
        sequences[row["copy_id"]] = seq

    checked_pairs = set()
    for i, row_a in enumerate(group):
        for row_b in group[i + 1:]:
            pair_key = tuple(sorted([row_a["copy_id"], row_b["copy_id"]]))
            if pair_key in checked_pairs:
                continue
            checked_pairs.add(pair_key)

            match_type = sequences_match(sequences[row_a["copy_id"]], sequences[row_b["copy_id"]])
            verdict = "assembly_duplicate" if match_type else "independent_copy"

            if match_type:
                # Prefer chromosome-style accession (e.g. CM/CP prefixes) over WGS contig accessions
                def is_chromosome_style(sseqid):
                    acc = clean_seqid(sseqid)
                    return acc.upper().startswith(("CM", "CP", "NC_"))

                a_is_chr = is_chromosome_style(row_a["sseqid"])
                b_is_chr = is_chromosome_style(row_b["sseqid"])
                if a_is_chr and not b_is_chr:
                    excluded_copy_ids.add(row_b["copy_id"])
                elif b_is_chr and not a_is_chr:
                    excluded_copy_ids.add(row_a["copy_id"])
                else:
                    # can't distinguish by naming convention -- flag for manual call, exclude neither automatically
                    verdict = "assembly_duplicate_ambiguous_which_to_keep"

            log_rows.append({
                "gene": gene, "accession": accession,
                "copy_id_a": row_a["copy_id"], "copy_id_b": row_b["copy_id"],
                "match_type": match_type or "no_match",
                "verdict": verdict,
                "excluded": row_b["copy_id"] if row_a["copy_id"] not in excluded_copy_ids and row_b["copy_id"] in excluded_copy_ids
                            else (row_a["copy_id"] if row_a["copy_id"] in excluded_copy_ids else ""),
            })

deduped_rows = [row for row in rows if row["copy_id"] not in excluded_copy_ids]

with open(OUT_FILE, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(deduped_rows)

log_fieldnames = ["gene", "accession", "copy_id_a", "copy_id_b", "match_type", "verdict", "excluded"]
with open(LOG_FILE, "w") as fh:
    writer = csv.DictWriter(fh, fieldnames=log_fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(log_rows)

print(f"Multi-locus (gene, genome) pairs checked: {len(log_rows)}")
print(f"Assembly duplicates auto-resolved: {sum(1 for r in log_rows if r['verdict'] == 'assembly_duplicate')}")
print(f"Ambiguous duplicates needing manual call: {sum(1 for r in log_rows if r['verdict'] == 'assembly_duplicate_ambiguous_which_to_keep')}")
print(f"Confirmed independent copies: {sum(1 for r in log_rows if r['verdict'] == 'independent_copy')}")
print(f"Loci excluded as duplicates: {len(excluded_copy_ids)}")
print(f"Final deduped homolog count: {len(deduped_rows)} (was {len(rows)})")
