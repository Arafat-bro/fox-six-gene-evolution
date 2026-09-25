import glob, os, collections
from Bio import SeqIO
from Bio.Seq import Seq

MINCOV = 0.7      # fraction of the query protein that must be covered
MINFRAC = 0.9     # gene must be present in this fraction of genomes
qlen = {r.id: len(r.seq) for r in SeqIO.parse("results/core/queries.faa", "fasta")}
os.makedirs("results/core/prot", exist_ok=True)
tsv = open("results/core/hits.tsv", "w")
tsv.write("genome\tog\tcontig\tstart\tend\tstrand\tidentity\tcoverage\n")

for gff in sorted(glob.glob("results/core/gff/*.gff")):
    n = os.path.basename(gff)[:-4]
    genome = SeqIO.to_dict(SeqIO.parse(f"data/genomes/raw/flat/{n}.fasta", "fasta"))
    mrna, cds = {}, collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "mRNA":
            t = a["Target"].split()
            mrna[a["ID"]] = dict(q=t[0], contig=f[0], s=int(f[3]), e=int(f[4]),
                                 strand=f[6], score=float(f[5]),
                                 ident=float(a["Identity"]),
                                 cov=(int(t[2]) - int(t[1]) + 1) / qlen[t[0]])
        elif f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]), int(f[4])))
    best = {}
    for mid, m in mrna.items():
        if m["q"] not in best or m["score"] > best[m["q"]][1]["score"]:
            best[m["q"]] = (mid, m)
    with open(f"results/core/prot/{n}.faa", "w") as out:
        for q, (mid, m) in best.items():
            if m["cov"] < MINCOV:
                continue
            seq = Seq("".join(str(genome[m["contig"]].seq[s - 1:e])
                              for s, e in sorted(cds[mid])))
            if m["strand"] == "-":
                seq = seq.reverse_complement()
            seq = seq[:len(seq) // 3 * 3]
            p = str(seq.translate()).rstrip("*").replace("*", "X")
            out.write(f">{n}\n{p}\n" if False else f">{q}\n{p}\n")
            tsv.write(f"{n}\t{q}\t{m['contig']}\t{m['s']}\t{m['e']}\t{m['strand']}\t{m['ident']}\t{m['cov']:.3f}\n")
tsv.close()

# group by gene, keep genes present in at least MINFRAC of genomes
genomes = sorted(os.path.basename(f)[:-4] for f in glob.glob("results/core/prot/*.faa"))
bygene = collections.defaultdict(dict)
for g in genomes:
    for r in SeqIO.parse(f"results/core/prot/{g}.faa", "fasta"):
        bygene[r.id][g] = str(r.seq)
os.makedirs("results/core/og", exist_ok=True)
kept = 0
for og, d in bygene.items():
    if len(d) >= MINFRAC * len(genomes):
        kept += 1
        with open(f"results/core/og/{og}.faa", "w") as o:
            for g, s in d.items():
                o.write(f">{g}\n{s}\n")
print(len(genomes), "genomes;", kept, "genes kept out of", len(bygene))