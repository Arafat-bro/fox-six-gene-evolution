#!/usr/bin/env python3
"""
Extract reference protein and coding-sequence (CDS) nucleotide queries for the
six SIX effector genes from their NCBI GenBank flatfiles.

Protein sequences are taken directly from each CDS feature's /translation
qualifier (the submitter's own annotation), rather than re-derived by
translating raw nucleotide sequence, so that codon_start and intron handling
match the original GenBank record exactly.

Input:  data/references/raw/<accession>.gb
Output: data/references/processed/<GENE>_protein.fasta
        data/references/processed/<GENE>_cds_nt.fasta
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition

RAW_DIR = Path("data/references/raw")
OUT_DIR = Path("data/references/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENE_MAP = {
    "MK906592.1": "SIX1",
    "MK906598.1": "SIX3",
    "MK906607.1": "SIX5",
    "GQ268954.1": "SIX7",
    "MK906667.1": "SIX10",
    "MW160867.1": "SIX12",
}

summary_rows = []

for accession, gene in GENE_MAP.items():
    gb_path = RAW_DIR / f"{accession}.gb"
    if not gb_path.exists():
        raise FileNotFoundError(f"Expected input file not found: {gb_path}")

    record = SeqIO.read(gb_path, "genbank")

    cds_features = [f for f in record.features if f.type == "CDS"]
    if len(cds_features) != 1:
        raise ValueError(f"{accession} ({gene}): expected exactly 1 CDS feature, found {len(cds_features)}")
    cds = cds_features[0]

    if "translation" not in cds.qualifiers:
        raise ValueError(f"{accession} ({gene}): CDS has no /translation qualifier")
    protein_seq = cds.qualifiers["translation"][0]

    cds_nt_seq = cds.location.extract(record.seq)

    n_exons = len(cds.location.parts)
    codon_start = cds.qualifiers.get("codon_start", ["1"])[0]
    is_partial = isinstance(cds.location.start, BeforePosition) or isinstance(cds.location.end, AfterPosition)

    protein_path = OUT_DIR / f"{gene}_protein.fasta"
    nt_path = OUT_DIR / f"{gene}_cds_nt.fasta"

    with open(protein_path, "w") as fh:
        fh.write(f">{gene}_{accession}_protein\n{protein_seq}\n")
    with open(nt_path, "w") as fh:
        fh.write(f">{gene}_{accession}_cds_nt\n{str(cds_nt_seq)}\n")

    summary_rows.append({
        "gene": gene, "accession": accession, "n_exons": n_exons,
        "codon_start": codon_start, "partial": is_partial,
        "aa_len": len(protein_seq), "nt_len": len(cds_nt_seq),
        "nt_div3": len(cds_nt_seq) % 3 == 0,
        "starts_M": protein_seq.startswith("M"),
    })

hdr = f"{'gene':6} {'accession':13} {'exons':6} {'codon_start':11} {'partial':8} {'aa_len':7} {'nt_len':7} {'nt%3==0':8} {'starts_M':9}"
print(hdr)
for r in summary_rows:
    print(f"{r['gene']:6} {r['accession']:13} {r['n_exons']:<6} {r['codon_start']:<11} "
          f"{str(r['partial']):8} {r['aa_len']:<7} {r['nt_len']:<7} {str(r['nt_div3']):8} {str(r['starts_M']):9}")
