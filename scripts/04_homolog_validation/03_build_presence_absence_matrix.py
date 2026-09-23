"""
Build the SIX presence/absence + copy-number matrix, plus a companion
pseudogene-annotation matrix, from homolog_classification_final.tsv.

Two outputs, deliberately kept separate rather than merged into one table:

  six_presence_absence_matrix.tsv
      rows = all genomes in the project's genome list (not just genomes with
      hits -- true absence is itself a data point and must appear as a row
      of zeros, not be silently omitted).
      columns = SIX1, SIX3, SIX5, SIX7, SIX10, SIX12
      cell = count of loci for that (genome, gene) with final_status ==
      "retained". This is the official copy-number matrix.

  six_pseudogene_annotation_matrix.tsv
      same shape, cell = count of loci with final_status ==
      "excluded_pseudogene" for that (genome, gene).
      This exists so a "0" in the main matrix can be checked against this
      table to distinguish TRUE absence (0 in both) from
      PSEUDOGENIZED-ONLY presence (0 retained, but >=1 here) -- collapsing
      these two states into a single "0" would erase exactly the
      complete/truncated/pseudogenized distinction the assignment requires.

Also prints any mismatches between the genome accession list and the
accessions actually appearing in the classification table, so a silently
missing or extra genome row is caught here rather than discovered later.

Input:  results/six_homologs/classified/homolog_classification_final.tsv
        scripts/genome_accessions.txt  (one accession per line; the
            project's full genome list, independent of hit results)
Output: results/six_presence_absence_matrix.tsv
        results/six_pseudogene_annotation_matrix.tsv
"""
import csv
from pathlib import Path
from collections import defaultdict

IN_FILE = Path("results/six_homologs/classified/homolog_classification_final.tsv")
ACCESSIONS_FILE = Path("scripts/02_genome_selection/final_accessions.txt")
PA_OUT = Path("results/six_presence_absence_matrix.tsv")
PSEUDO_OUT = Path("results/six_pseudogene_annotation_matrix.tsv")

GENES = ["SIX1", "SIX3", "SIX5", "SIX7", "SIX10", "SIX12"]

with open(ACCESSIONS_FILE) as fh:
    genome_accessions = [line.strip() for line in fh if line.strip()]

with open(IN_FILE) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

accessions_in_table = {r["accession"] for r in rows}
missing_from_table = set(genome_accessions) - accessions_in_table
extra_in_table = accessions_in_table - set(genome_accessions)

if missing_from_table:
    print(f"WARNING: {len(missing_from_table)} genome(s) in accession list have "
          f"no rows at all in the classification table (i.e. zero hits across "
          f"all six genes for that genome -- confirm this is a true negative, "
          f"not a missing/failed BLAST run):")
    for acc in sorted(missing_from_table):
        print(f"  {acc}")

if extra_in_table:
    print(f"WARNING: {len(extra_in_table)} accession(s) appear in the "
          f"classification table but not in the genome accession list "
          f"(possible typo or stale entry):")
    for acc in sorted(extra_in_table):
        print(f"  {acc}")

retained_counts = defaultdict(int)
pseudogene_counts = defaultdict(int)

for row in rows:
    key = (row["accession"], row["gene"])
    if row["final_status"] == "retained":
        retained_counts[key] += 1
    elif row["final_status"] == "excluded_pseudogene":
        pseudogene_counts[key] += 1

def write_matrix(path, counts):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome"] + GENES)
        for acc in genome_accessions:
            writer.writerow([acc] + [counts.get((acc, gene), 0) for gene in GENES])

write_matrix(PA_OUT, retained_counts)
write_matrix(PSEUDO_OUT, pseudogene_counts)

print()
print(f"Presence/absence matrix written: {PA_OUT} ({len(genome_accessions)} genomes x {len(GENES)} genes)")
print(f"Pseudogene annotation matrix written: {PSEUDO_OUT}")

print()
print("Per-gene totals (retained copies, summed across all genomes):")
for gene in GENES:
    total = sum(retained_counts.get((acc, gene), 0) for acc in genome_accessions)
    genomes_with_gene = sum(1 for acc in genome_accessions if retained_counts.get((acc, gene), 0) > 0)
    print(f"  {gene}: {total} total copies, present in {genomes_with_gene}/{len(genome_accessions)} genomes")
