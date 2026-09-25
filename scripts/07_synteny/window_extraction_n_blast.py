mkdir -p results/synteny/windows

# 1. Build window coordinates from the linked pairs only
python3 << 'EOF'
import pandas as pd
df = pd.read_csv("results/synteny/six7_10_12_linkage_table.tsv", sep="\t")
linked = df[df["same_contig_as_partner"] == "yes"]

with open("results/synteny/windows/window_coords.tsv", "w") as f:
    for (genome, contig), grp in linked.groupby(["genome", "contig"]):
        start = max(0, grp["start"].min() - 5000)
        end = grp["end"].max() + 5000
        genes = "_".join(sorted(grp["gene"].str.replace(r"_copy\d+", "", regex=True).unique()))
        f.write(f"{genome}\t{contig}\t{start}\t{end}\t{genes}\n")
        print(genome, contig, start, end, genes)
EOF

# 2. Confirm your genome FASTA naming before extracting (quick sanity check)
ls data/genomes/raw/flat/ | head -3

# 3. Extract each window
while IFS=$'\t' read -r genome contig start end genes; do
  fasta="data/genomes/raw/flat/${genome}.fasta"
  [ -f "$fasta" ] || { echo "MISSING FASTA for $genome"; continue; }
  seqkit faidx "$fasta" "${contig}:${start}-${end}" > "results/synteny/windows/${genome}_${genes}.fasta"
done < results/synteny/windows/window_coords.tsv

ls -la results/synteny/windows/

# 4. BLAST the Fol4287 reference window against every other window that shares the same gene pair
ref=$(ls results/synteny/windows/GCA_003315725.1_*.fasta 2>/dev/null)
echo "Reference window: $ref"
for f in results/synteny/windows/*.fasta; do
  [ "$f" = "$ref" ] && continue
  echo "=== Fol4287 window vs $(basename $f) ==="
  blastn -query "$ref" -subject "$f" -outfmt "6 pident length mismatch gapopen evalue bitscore"
  echo
done