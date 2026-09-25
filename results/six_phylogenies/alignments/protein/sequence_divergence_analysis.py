from Bio import AlignIO
from itertools import combinations
import pandas as pd
import os
BASE = "/home/arafat/Foxy_SIX_Evolution/results/six_phylogenies/alignments/protein"
files = {
    "SIX1": f"{BASE}/SIX1.aln.fasta",
    "SIX3": f"{BASE}/SIX3.aln.fasta",
    "SIX5": f"{BASE}/SIX5.aln.fasta",
    "SIX7": f"{BASE}/SIX7.aln.fasta",
    "SIX10": f"{BASE}/SIX10.aln.fasta",
    "SIX12": f"{BASE}/six12_aligned.fasta"
}


# Put your aligned FASTA files here
files = {
    "SIX1": "SIX1.aln.fasta",
    "SIX3": "SIX3.aln.fasta",
    "SIX5": "SIX5.aln.fasta",
    "SIX7": "SIX7.aln.fasta",
    "SIX10": "SIX10.aln.fasta",
    "SIX12": "six12_aligned.fasta"
}


def calculate_identity(seq1, seq2):
    """
    Calculate pairwise amino acid identity
    ignoring gaps
    """
    matches = 0
    compared = 0

    for a, b in zip(seq1, seq2):

        if a == "-" or b == "-":
            continue

        compared += 1

        if a == b:
            matches += 1

    if compared == 0:
        return 0

    return (matches / compared) * 100



results = []


for gene, filename in files.items():

    print("Processing:", gene)

    alignment = AlignIO.read(filename, "fasta")

    sequences = [str(record.seq) for record in alignment]

    n_sequences = len(sequences)

    alignment_length = alignment.get_alignment_length()


    # Count conserved and variable sites
    conserved = 0
    variable = 0


    for column in range(alignment_length):

        site = [
            seq[column]
            for seq in sequences
            if seq[column] != "-"
        ]

        if len(set(site)) <= 1:
            conserved += 1
        else:
            variable += 1



    # Pairwise identity
    identities = []

    for s1, s2 in combinations(sequences, 2):

        identity = calculate_identity(s1, s2)

        identities.append(identity)



    mean_identity = sum(identities) / len(identities)

    divergence = 100 - mean_identity



    results.append({

        "Gene": gene,
        "Number_of_sequences": n_sequences,
        "Alignment_length_aa": alignment_length,
        "Conserved_sites": conserved,
        "Variable_sites": variable,
        "Mean_pairwise_identity_%": round(mean_identity,2),
        "Mean_divergence_%": round(divergence,2)

    })



df = pd.DataFrame(results)


print("\nSequence Divergence Summary\n")

print(df)


df.to_csv(
    "SIX_sequence_divergence_summary.csv",
    index=False
)

print("\nSaved: SIX_sequence_divergence_summary.csv")