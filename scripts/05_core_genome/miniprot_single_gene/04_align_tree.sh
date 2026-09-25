#!/bin/bash
mkdir -p results/core/aln results/core/trim
for f in results/core/og/*.faa; do
  og=$(basename $f .faa)
  [ -s results/core/trim/$og.faa ] && continue
  mafft --auto --thread 8 $f > results/core/aln/$og.faa 2>/dev/null
  trimal -in results/core/aln/$og.faa -out results/core/trim/$og.faa -automated1
done

python - <<'EOF'
import glob, os
from Bio import SeqIO
genomes = sorted(os.path.basename(f)[:-4] for f in glob.glob("results/core/prot/*.faa"))
seqs = {g: [] for g in genomes}
for f in sorted(glob.glob("results/core/trim/*.faa")):
    aln = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    L = len(next(iter(aln.values())))
    if L < 50:
        continue
    for g in genomes:
        seqs[g].append(aln.get(g, "-" * L))
with open("results/core/concat.faa", "w") as o:
    for g in genomes:
        o.write(f">{g}\n{''.join(seqs[g])}\n")
print("concatenated length:", len("".join(seqs[genomes[0]])))
EOF

which iqtree3 iqtree
iqtree3 -s results/core/concat.faa -m MFP -B 1000 -T 8 -mem 6G --prefix results/core/core_tree