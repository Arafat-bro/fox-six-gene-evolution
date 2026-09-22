"""
For every locus in genuine_homologs_deduped_final.tsv AND
borderline_for_manual_review.tsv:
  1. Extract the genomic sequence at its coordinates, oriented to the coding
     strand (reverse-complemented if sstrand == minus).
  2. Translate the extracted sequence in all three forward reading frames of
     that oriented sequence (frame selection is necessary because merged-HSP
     boundaries from a protein-vs-translated-genome search are not guaranteed
     to fall on a codon boundary).
  3. Select the frame with the longest stop-free ORF fragment; ties broken by
     closeness of that ORF's length to the reference protein's length for the
     same gene.
  4. Classify the locus as complete / truncated / pseudogenized based on the
     winning frame's internal-stop count and length relative to the reference
     protein.
  5. Independently flag loci whose coordinates sit at or near the scaffold's
     start/end as assembly-boundary-affected -- a caveat layered on top of,
     not a replacement for, the translation-based classification.

Both input tables are processed in the same pass and tagged with a
"source_table" column ("genuine" or "borderline") so borderline loci remain
visibly distinct in the output rather than being silently merged into the
main homolog count. Borderline loci get the same translation evidence so
their manual accept/reject calls can be made on ORF quality, not identity
percentage alone.

Input:  results/six_homologs/filtered/genuine_homologs_deduped_final.tsv
        results/six_homologs/filtered/borderline_for_manual_review.tsv
        data/genomes/raw/flat/<accession>.fasta
        data/references/processed/<gene>_protein.fasta
Output: results/six_homologs/classified/homolog_classification.tsv
        results/six_homologs/classified/translated_proteins/<gene>/<copy_id>.faa
        results/six_homologs/classified/nt_sequences/<gene>/<copy_id>.fasta
"""
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GENUINE_FILE = Path("results/six_homologs/filtered/genuine_homologs_deduped_final.tsv")
BORDERLINE_FILE = Path("results/six_homologs/filtered/borderline_for_manual_review.tsv")
GENOME_DIR = Path("data/genomes/raw/flat")
REF_PROTEIN_DIR = Path("data/references/processed")
OUT_FILE = Path("results/six_homologs/classified/homolog_classification.tsv")
PROT_OUT_DIR = Path("results/six_homologs/classified/translated_proteins")
NT_OUT_DIR = Path("results/six_homologs/classified/nt_sequences")

BOUNDARY_MARGIN_BP = 10          # locus edge within this many bp of scaffold edge -> boundary flag
TRUNCATION_THRESHOLD = 0.90      # protein_length / ref_length below this -> "truncated" instead of "complete"


def clean_seqid(sseqid):
    parts = sseqid.split("|")
    return parts[1] if len(parts) > 1 else sseqid


def load_ref_lengths():
    lengths = {}
    for fasta in REF_PROTEIN_DIR.glob("*_protein.fasta"):
        gene = fasta.stem.replace("_protein", "")
        rec = next(SeqIO.parse(fasta, "fasta"))
        lengths[gene] = len(str(rec.seq).rstrip("*"))
    return lengths


def best_frame(nt_seq):
    results = []
    for frame in range(3):
        sub = nt_seq[frame:]
        sub = sub[: len(sub) - (len(sub) % 3)]
        translation = str(Seq(sub).translate())
        fragments = translation.split("*")
        longest = max(fragments, key=len)
        results.append((frame, translation, len(longest)))
    return results


def classify(translation, ref_len):
    body = translation.rstrip("*")
    internal_stops = body.count("*")
    ends_in_stop = translation.endswith("*")
    protein_len = len(body)
    if internal_stops > 0:
        return "pseudogenized", internal_stops, ends_in_stop, protein_len
    if ref_len and protein_len < TRUNCATION_THRESHOLD * ref_len:
        return "truncated", internal_stops, ends_in_stop, protein_len
    return "complete", internal_stops, ends_in_stop, protein_len


def load_rows(path, source_label):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["source_table"] = source_label
    return rows


ref_lengths = load_ref_lengths()
rows = load_rows(GENUINE_FILE, "genuine") + load_rows(BORDERLINE_FILE, "borderline")

genome_cache = {}
out_rows = []

for row in rows:
    accession = row["accession"]
    gene = row["gene"]
    copy_id = row["copy_id"]

    if accession not in genome_cache:
        fasta_path = GENOME_DIR / f"{accession}.fasta"
        genome_cache[accession] = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    genome_records = genome_cache[accession]

    scaffold_acc = clean_seqid(row["sseqid"])
    scaffold_rec = genome_records.get(scaffold_acc)
    if scaffold_rec is None:
        out_rows.append({**row, "error": f"scaffold {scaffold_acc} not found in genome fasta"})
        continue
    scaffold_len = len(scaffold_rec.seq)

    start, end = int(row["genomic_start"]), int(row["genomic_end"])
    lo, hi = min(start, end), max(start, end)
    raw_seq = str(scaffold_rec.seq[lo - 1:hi])

    if row["sstrand"].lower().startswith("minus"):
        oriented_seq = str(Seq(raw_seq).reverse_complement())
    else:
        oriented_seq = raw_seq

    ref_len = ref_lengths.get(gene)
    frame_results = best_frame(oriented_seq)

    def score(fr):
        frame, translation, orf_len = fr
        dist = abs(orf_len - ref_len) if ref_len else 0
        return (-orf_len, dist)

    frame, translation, orf_len = sorted(frame_results, key=score)[0]

    classification, internal_stops, ends_in_stop, protein_len = classify(translation, ref_len)

    boundary_flag = (lo <= BOUNDARY_MARGIN_BP) or (scaffold_len - hi <= BOUNDARY_MARGIN_BP)

    gene_prot_dir = PROT_OUT_DIR / gene
    gene_nt_dir = NT_OUT_DIR / gene
    gene_prot_dir.mkdir(parents=True, exist_ok=True)
    gene_nt_dir.mkdir(parents=True, exist_ok=True)
    with open(gene_prot_dir / f"{copy_id}.faa", "w") as fh:
        fh.write(f">{copy_id} frame={frame} source={row['source_table']}\n{translation.rstrip('*')}\n")
    with open(gene_nt_dir / f"{copy_id}.fasta", "w") as fh:
        fh.write(f">{copy_id}\n{oriented_seq}\n")

    out_rows.append({
        **row,
        "scaffold_length": scaffold_len,
        "frame_selected": frame,
        "protein_length": protein_len,
        "ref_protein_length": ref_len,
        "internal_stop_count": internal_stops,
        "ends_in_stop": ends_in_stop,
        "assembly_boundary_flag": boundary_flag,
        "classification": classification,
    })

fieldnames = list(rows[0].keys()) + [
    "scaffold_length", "frame_selected", "protein_length", "ref_protein_length",
    "internal_stop_count", "ends_in_stop", "assembly_boundary_flag", "classification",
]
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_FILE, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Loci processed: {len(out_rows)}")
for src in ["genuine", "borderline"]:
    print(f"  source={src}: {sum(1 for r in out_rows if r.get('source_table') == src)}")
for cat in ["complete", "truncated", "pseudogenized"]:
    print(f"  {cat}: {sum(1 for r in out_rows if r.get('classification') == cat)}")
print(f"  assembly_boundary_flag=True: {sum(1 for r in out_rows if r.get('assembly_boundary_flag'))}")
