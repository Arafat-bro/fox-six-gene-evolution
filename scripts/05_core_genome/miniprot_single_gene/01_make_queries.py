import glob, os
from Bio import SeqIO

BUSCO = "results/busco"
REF = "GCA_003615085.1"   # genome whose proteins are used as queries
SUB = "run_hypocreales_odb12.2/busco_sequences/single_copy_busco_sequences"

sc = {}
for d in glob.glob(f"{BUSCO}/*/{SUB}"):
    g = d.split("/")[2]
    sc[g] = {os.path.basename(f)[:-4] for f in glob.glob(d + "/*.faa")}
print(len(sc), "BUSCO folders found:", sorted(sc))

common = set.intersection(*sc.values())
print(len(common), "single copy genes shared by all", len(sc), "genomes")

os.makedirs("results/core", exist_ok=True)
with open("results/core/queries.faa", "w") as out:
    for og in sorted(common):
        rec = next(SeqIO.parse(f"{BUSCO}/{REF}/{SUB}/{og}.faa", "fasta"))
        out.write(f">{og}\n{str(rec.seq).replace('*', '')}\n")