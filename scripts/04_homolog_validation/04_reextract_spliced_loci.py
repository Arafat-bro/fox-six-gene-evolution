"""
Re-extract and translate the 5 loci belonging to genuinely spliced reference
genes (SIX5, SIX10) where the naive min-to-max genomic span used by
03_merge_hsps_into_loci.py would incorrectly include intron sequence,
producing a frameshifted, garbage translation.

For each affected locus: extract a padded genomic window around the naive
span, run exonerate (protein2genome model) with the reference protein as
query, parse the predicted exon coordinates from its GFF output, splice
together only the exon sequence (skipping introns), then translate and
classify that corrected sequence exactly as 01_extract_translate_classify.py
does for ordinary loci.

Reference confirmed spliced via GenBank CDS join() records:
  SIX5  (MK906607.1): join(57..100,166..275,328..395,445..582)
  SIX10 (MK906667.1): join(21..100,171..540)
SIX1 and SIX7 confirmed NOT spliced (plain CDS ranges, no join()) and are
therefore excluded from this reprocessing -- their existing classification
results from the naive span extraction stand as-is.

Input:  results/six_homologs/classified/homolog_classification_final.tsv
        data/genomes/raw/flat/<accession>.fasta
        data/references/processed/<gene>_protein.fasta
Output: results/six_homologs/classified/homolog_classification_final.tsv
        (patched in place for just these 5 rows)
        results/six_homologs/classified/translated_proteins/<gene>/<copy_id>.faa
        (overwritten for these 5 loci)
        docs/decisions/spliced_locus_reextraction_log.tsv
"""
import csv
import subprocess
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GENOME_DIR = Path("data/genomes/raw/flat")
REF_PROTEIN_DIR = Path("data/references/processed")
CLASS_FILE = Path("results/six_homologs/classified/homolog_classification_final.tsv")
PROT_OUT_DIR = Path("results/six_homologs/classified/translated_proteins")
LOG_FILE = Path("docs/decisions/spliced_locus_reextraction_log.tsv")
WORKDIR = Path("results/six_homologs/classified/exonerate_tmp")
WORKDIR.mkdir(parents=True, exist_ok=True)

PADDING_BP = 1000
TRUNCATION_THRESHOLD = 0.90

AFFECTED_COPY_IDS = {
    "SIX5_GCA_003315725.1_copy1",
    "SIX5_GCA_003977725.1_copy1",
    "SIX10_GCA_003315725.1_copy1",
    "SIX10_GCA_003977725.1_copy1",
    "SIX10_GCA_013423245.1_copy1",
}


def clean_seqid(sseqid):
    parts = sseqid.split("|")
    return parts[1] if len(parts) > 1 else sseqid


def run_exonerate(query_fasta, target_fasta):
    cmd = [
        "exonerate", "--model", "protein2genome",
        "--query", str(query_fasta), "--target", str(target_fasta),
        "--showtargetgff", "yes", "--showalignment", "no", "--showvulgar", "no",
        "--bestn", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout


def parse_exon_coords(gff_text):
    exons = []
    for line in gff_text.splitlines():
        if "\texon\t" in line:
            fields = line.split("\t")
            start, end = int(fields[3]), int(fields[4])
            exons.append((start, end))
    return sorted(exons)


with open(CLASS_FILE) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
    fieldnames = list(rows[0].keys())

log_rows = []

for row in rows:
    if row["copy_id"] not in AFFECTED_COPY_IDS:
        continue

    gene, accession = row["gene"], row["accession"]
    genome_records = SeqIO.to_dict(SeqIO.parse(GENOME_DIR / f"{accession}.fasta", "fasta"))
    scaffold_acc = clean_seqid(row["sseqid"])
    scaffold_rec = genome_records[scaffold_acc]
    scaffold_len = len(scaffold_rec.seq)

    naive_lo, naive_hi = int(row["genomic_start"]), int(row["genomic_end"])
    window_lo = max(1, naive_lo - PADDING_BP)
    window_hi = min(scaffold_len, naive_hi + PADDING_BP)
    window_seq = str(scaffold_rec.seq[window_lo - 1:window_hi])

    target_fasta = WORKDIR / f"{row['copy_id']}_window.fasta"
    with open(target_fasta, "w") as fh:
        fh.write(f">{row['copy_id']}_window\n{window_seq}\n")

    query_fasta = REF_PROTEIN_DIR / f"{gene}_protein.fasta"
    gff_text = run_exonerate(query_fasta, target_fasta)
    exon_coords = parse_exon_coords(gff_text)

    if not exon_coords:
        log_rows.append({"copy_id": row["copy_id"], "status": "exonerate_found_no_exons"})
        continue

    exon_seq = "".join(window_seq[s - 1:e] for s, e in exon_coords)

    # exonerate's GFF exon coordinates are on the strand it aligned to; if the
    # sequence came back needing reverse-complementing (i.e. the CDS reads
    # antisense relative to how we extracted the window), detect via presence
    # of a valid ORF in either orientation and pick the one with fewer stops.
    def translate_and_score(nt_seq):
        best = None
        for frame in range(3):
            sub = nt_seq[frame:]
            sub = sub[: len(sub) - (len(sub) % 3)]
            t = str(Seq(sub).translate())
            stops = t.rstrip("*").count("*")
            if best is None or stops < best[2]:
                best = (frame, t, stops)
        return best

    fwd_frame, fwd_t, fwd_stops = translate_and_score(exon_seq)
    rev_frame, rev_t, rev_stops = translate_and_score(str(Seq(exon_seq).reverse_complement()))
    if rev_stops < fwd_stops:
        frame, translation, internal_stops = rev_frame, rev_t, rev_stops
        orientation = "reverse_complement"
    else:
        frame, translation, internal_stops = fwd_frame, fwd_t, fwd_stops
        orientation = "forward"

    ref_len = len(str(next(SeqIO.parse(query_fasta, "fasta")).seq).rstrip("*"))
    body = translation.rstrip("*")
    protein_len = len(body)
    if internal_stops > 0:
        classification = "pseudogenized"
    elif protein_len < TRUNCATION_THRESHOLD * ref_len:
        classification = "truncated"
    else:
        classification = "complete"

    row["frame_selected"] = f"exonerate:{orientation}:{frame}"
    row["protein_length"] = protein_len
    row["internal_stop_count"] = internal_stops
    row["ends_in_stop"] = translation.endswith("*")
    row["classification"] = classification
    row["final_status"] = "excluded_pseudogene" if classification == "pseudogenized" else "retained"

    gene_prot_dir = PROT_OUT_DIR / gene
    with open(gene_prot_dir / f"{row['copy_id']}.faa", "w") as fh:
        fh.write(f">{row['copy_id']} exonerate_reextracted n_exons={len(exon_coords)}\n{body}\n")

    log_rows.append({
        "copy_id": row["copy_id"], "gene": gene, "accession": accession,
        "n_exons_found": len(exon_coords), "exon_coords_in_window": str(exon_coords),
        "orientation": orientation, "protein_length": protein_len,
        "internal_stop_count": internal_stops, "classification": classification,
        "status": "reextracted_ok",
    })

with open(CLASS_FILE, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

with open(LOG_FILE, "w", newline="") as fh:
    log_fields = ["copy_id", "gene", "accession", "n_exons_found", "exon_coords_in_window",
                  "orientation", "protein_length", "internal_stop_count", "classification", "status"]
    writer = csv.DictWriter(fh, fieldnames=log_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for r in log_rows:
        writer.writerow(r)

print(f"Re-extracted {len(log_rows)} loci.")
for r in log_rows:
    print(f"  {r['copy_id']}: {r.get('classification', r['status'])} "
          f"(exons={r.get('n_exons_found','?')}, stops={r.get('internal_stop_count','?')})")
