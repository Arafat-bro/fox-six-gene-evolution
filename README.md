# Foxy_SIX_Evolution — SIX Effector Gene Evolution in *Fusarium oxysporum*

Reconstructing the evolutionary history of six SIX (Secreted In Xylem)
effector genes — SIX1, SIX3, SIX5, SIX7, SIX10, SIX12 — across 27
*Fusarium oxysporum* genome assemblies, and testing whether each gene's
pattern is better explained by vertical inheritance, duplication/loss,
sequence divergence, or movement of an accessory genomic region.

Repo: https://github.com/Arafat-bro/fox-six-gene-evolution

---

## 1. The argument this pipeline builds toward

1. A **core-genome tree**, built from conserved single-copy orthologs, is
   the organismal reference — "what actually happened to these genomes."
2. Each SIX gene gets its **own tree** from its retained homologs.
3. Each SIX tree is **compared to the core tree**. Agreement → consistent
   with vertical inheritance. Disagreement → duplication/loss, divergence,
   or possible horizontal movement — but incongruence alone never proves
   transfer.
4. The **SIX7–SIX10–SIX12 synteny check** is the corroborating evidence:
   if the three genes stay linked, in order, with conserved flanking
   sequence, across genomes that are *not* close on the core tree, that's
   much stronger evidence for block-level movement than any single gene's
   tree position.

Everything below produces one piece of that argument. Result numbers here
are taken directly from the final report (`res_updated.docx`) — where the
report and this README ever disagree, the report is the source of truth.

---

## 2. Pipeline, in run order

### Step 1 — Genome selection + metadata
27 *F. oxysporum* assemblies selected from NCBI GenBank and curated into
`results/genome_metadata.tsv` (isolate, forma specialis, host, geographic
origin, accession, assembly level, genome size, contig count, N50,
sequencing technology).

**Why these 27 and not a mechanical "top N" pull:** the panel had to
support the comparisons the biological questions actually ask —
multiple independent host/forma-specialis lineages, a mix of assembly
qualities acknowledged explicitly (not excluded, since SIX genes sit on
accessory chromosomes that fragment easily and absence calls are only
trustworthy once assembly quality is accounted for), the nonpathogenic
isolate **170** and pathogenicity-chromosome reference isolates
(**Fol4287**, **race3**) needed as anchors, and enough isolates per forma
specialis to distinguish "fixed within this lineage" from "only sampled
once." Host/forma-specialis/geography metadata is not present in NCBI's
summary tables in a clean form, so this was curated by hand from each
genome's BioSample record — logged as a manual step, not scripted.

### Step 2 — Core-genome marker discovery
`scripts/05_core_genome/miniprot_single_gene/01_make_queries.py`
BUSCO (`hypocreales_odb12.2`) was run on 6 representative genomes to
identify a shared single-copy candidate marker set.

**Why BUSCO on 6 genomes, not all 27:** running BUSCO — which does *ab
initio* gene prediction — on all 27 assemblies is computationally
expensive and unnecessary just to define *which* loci are usable markers.
Six representative genomes are enough to discover a robust shared
single-copy set; the more expensive step (checking whether each candidate
marker is present and single-copy in *every* genome) is handed off to a
much faster tool in the next step.

**Why `hypocreales_odb12.2` and not the broader `fungi_odb10`:**
Hypocreales is a much narrower, more taxonomically relevant lineage
dataset for *Fusarium* than the whole-fungal-kingdom set — a narrower
reference set gives more precisely conserved, single-copy markers for
this genus specifically, rather than markers that are merely conserved
across all fungi.

### Step 3 — Core-genome ortholog mapping
`scripts/05_core_genome/miniprot_single_gene/02_run_miniprot.sh`
The candidate marker proteins from Step 2 were mapped against all 27
genomes using **Miniprot v0.18** (spliced protein-to-genome alignment).

**Why Miniprot here instead of re-running BUSCO on all 27 genomes, or
using tblastn:** Miniprot is built specifically for aligning a known
protein set back onto a genome (including across introns), which is
exactly this problem — "does genome X have this exact marker, and where."
It is far cheaper than running full BUSCO gene prediction on every
genome, which matters on the hardware actually used (an 8 GB RAM, 8-core
workstation — not a cluster). tblastn was deliberately reserved for the
SIX-gene search (Step 6), where raw six-frame translated search against
unannotated sequence is the point; core-gene mapping doesn't need that,
since the marker proteins are already well-defined and mostly intron-poor
relative to effector loci.

### Step 4 — Ortholog filtering
`scripts/05_core_genome/miniprot_single_gene/03_extract.py`
Filtering thresholds: `MINCOV=0.7`, `MINFRAC=0.9`. Markers with an
internal stop codon in the predicted amino acid sequence were also
removed — these were identified as a recurring artifact of the gene
prediction step itself, not a biological signal, so keeping them would
have injected noise into the alignment rather than real divergence.
After filtering, **3,620 shared high-confidence orthologs** were retained
for the final tree (this is the number reported in the final analysis;
an earlier filtering pass had retained 3,911 of 3,912 candidates before
the internal-stop-codon cleanup was added).

### Step 5 — Core-genome alignment and tree
`scripts/05_core_genome/miniprot_single_gene/04_align_tree.sh`
MAFFT (`--auto`) alignment per marker → trimAl (`-automated1`) trimming →
concatenation into a supermatrix → IQ-TREE 3, best-fit model selected by
ModelFinder (LG+I+G), 1000 ultrafast bootstrap replicates.
Output: `results/core/core_tree_fast.treefile`.

**Why the tree itself was built from a defined subset of genomes (Fol4287,
Fo47, one representative forma-specialis isolate, and three additional
genomes) rather than all 27:** concatenated supermatrix alignment and
tree inference over the full panel was not feasible on the available
hardware within the project's time budget. The subset was chosen to keep
the isolates the rest of the analysis depends on (the reference strain,
the nonpathogenic control, and enough lineage spread for the tree to be
informative) rather than picking genomes at random.

**Why concatenation instead of a coalescent method (e.g. ASTRAL):**
concatenation is a defensible simplification for a project at this scope;
ASTRAL would more explicitly model gene-tree discordance from incomplete
lineage sorting, but was not pursued given the same time/compute
constraints. This trade-off is stated directly, not hidden.

### Step 6 — SIX gene homolog search
tblastn of the six reference SIX proteins (SIX1 = MK906592.1, SIX3 =
MK906598.1, SIX5 = MK906607.1, SIX7 = GQ268954.1, SIX10 = MK906667.1,
SIX12 = MW160867.1) against all 27 genome assemblies.
Thresholds: e-value ≤ 1e-10, query coverage ≥ 70%.
Outputs: `results/six_presence_absence_matrix.tsv`,
`results/six_pseudogene_annotation_matrix.tsv`.

**Why tblastn and not annotation lookup or blastn:** SIX genes are short,
fast-evolving, lineage-specific effectors that annotation pipelines
routinely miss. tblastn searches the raw six-frame-translated assembly
directly and tolerates codon-level divergence better than a nucleotide
(blastn) search, so it finds real homologs even where nothing was
annotated.

Each retained hit was manually classified as **complete / truncated /
pseudogenized** (premature stop codon or frameshift), and hits touching a
contig edge were flagged separately — an apparent absence there is
indistinguishable from an assembly gap without that flag. Borderline
hits (identity/coverage near the threshold) were flagged for manual
review rather than silently kept or dropped.

### Step 7 — Per-SIX-gene phylogenies
Protein-based trees were built for **SIX1, SIX3, SIX7, and SIX10** —
these had enough retained, non-identical homologs to support tree
inference (SIX1: 21 retained copies across 12 genomes; SIX3: 3; SIX7: 6
across 3 genomes; SIX10: 5 across 4 genomes).

**SIX5 and SIX12 were explicitly excluded from phylogenetic
reconstruction** — SIX5 has only 2 retained, essentially identical copies
(100% pairwise identity), and SIX12 has only 1 retained functional copy
(the rest are pseudogenized) — neither supports a meaningful tree. This
is a scoping decision, not an oversight, and is reported as such.

### Step 8 — SIX-tree vs. core-tree comparison
Comparison done with **DendroPy**, matching genomes by accession-based
labels and treating duplicate copies within a genome separately for
copy-number purposes.
- **SIX1**: enough shared taxa (12 genomes) for a formal **Robinson–Foulds
  distance** comparison against the core tree.
- **SIX3, SIX7, SIX10**: only 3 genomes each — too few for a meaningful RF
  distance, so comparison here is by topology/branch-length inspection
  only, and that limitation is stated explicitly rather than reporting a
  misleadingly precise statistic on 3 tips.

### Step 9 — SIX7–SIX10–SIX12 synteny analysis
For every genome carrying ≥2 of the three genes, physical linkage, gene
order, orientation, and intergenic distance were compared, and flanking
regions were checked with pairwise sequence alignment (not clinker — see
Limitations below) to test conservation of the whole block rather than
individual gene similarity alone.

**Key findings** (see report §3, Figure 6): Fol4287 carries the region as
SIX10–SIX12–SIX7 on one contig; race3 carries the same three genes in
reverse order and orientation with matching intergenic spacing, and the
two regions are 99.94% identical across the full ~15 kb window including
flanking sequence — the strongest single result in the project. The
*lini* isolate retains a linked SIX7–SIX12–SIX10 block with different
spacing, plus a second SIX7/SIX10 region ~3.7 Mb away lacking SIX12
(evidence of lineage-specific duplication/rearrangement). The *cepae*
isolate retains only a partial cluster (SIX12–SIX10 linked; no SIX7
copy detected there).

---

## 3. Environment

Local workstation, 8 GB RAM, 8 cores, WSL2 Ubuntu. Conda environment
`six_project`:

| Tool | Version |
|---|---|
| Python | 3.13.15 |
| Miniprot | 0.18 |
| BUSCO | 6.0.0 |
| MAFFT | 7.526 |
| trimAl | 1.5.1 |
| IQ-TREE | 3.1.3 |
| BLAST+, seqkit, Biopython, DendroPy | see `environment.yml` for pinned versions |

---

## 4. Manual (non-scripted) steps — logged explicitly

- Genome selection criteria and forma specialis/host/geography metadata:
  curated by hand from NCBI BioSample records
  (`results/genome_metadata.tsv`) — NCBI's summary tables don't carry this
  cleanly.
- Reference SIX gene GenBank records: CDS boundaries, strand, and
  mRNA-vs-genomic status manually checked before use as BLAST queries.
- Borderline BLAST hits: flagged for manual review rather than
  automatically included or excluded.
- SIX complete/truncated/pseudogenized classification: manually confirmed
  against premature stop codons and frameshifts.

---

## 5. Limitations (stated in the report, not omitted here)

- **dN/dS (Ka/Ks) analysis was scoped out.** The retained ortholog sets
  for most SIX genes were too small and, in several cases (e.g. SIX5's
  two copies at 100% identity), too nearly identical to yield an
  interpretable codon-based selection signal within the project's time
  constraints.
- **SIX3, SIX7, and SIX10 trees are 3-tip trees.** Bootstrap values there
  confirm the observed pairwise clustering is stable but cannot test
  broader alternative topologies — there simply aren't enough taxa.
- **SIX5 and SIX12 have no phylogenetic tree at all** — too few
  non-identical sequences to justify one (see Step 7).
- **The core-genome tree was built from a 6-genome subset**, not the full
  27-genome panel, due to compute/time constraints — this is a
  simplification, not a hidden limitation.
- **Gene-tree incongruence is not, by itself, evidence of horizontal
  transfer** — the report treats the Fol4287–race3 SIX7/SIX10/SIX12
  synteny result as the strongest corroborating evidence precisely
  because incongruence alone (seen at SIX1, SIX3, and the SIX12
  pseudogene as well) cannot distinguish transfer from duplication/loss
  or ancestral retention.
- Synteny comparison used pairwise alignment of flanking regions rather
  than a dedicated visualization tool like `clinker` — sufficient to
  establish conservation/order/spacing, but without the polished
  multi-genome comparative figure `clinker` would produce.

---

## 6. Results summary (see `res_updated.docx` for full detail)

| Gene | Genomes w/ retained copies | Retained copies | Complete | Truncated | Identity range (%) |
|---|---|---|---|---|---|
| SIX1 | 12 | 21 | 18 | 3 | 64.1–100.0 |
| SIX3 | 3 | 3 | 3 | 0 | 86.4–100.0 |
| SIX5 | 2 | 2 | 2 | 0 | 84.3 |
| SIX7 | 3 | 6 | 5 | 1 | 85.7–100.0 |
| SIX10 | 4 | 5 | 3 | 2 | 90.1–100.0 |
| SIX12 | 1 (functional) | 1 | 1 | 0 | 92.5 |

**Headline conclusion:** the six SIX genes do not follow one evolutionary
model. SIX1 shows broad distribution driven by duplication and
lineage-specific loss, largely independent of host association. SIX3,
SIX5, SIX7, SIX10, and SIX12 are restricted to three host-associated
lineages (tomato: Fol4287/race3; onion: FoC_Fus2; flax: *lini*) whose
core-genome positions do not form a clade — meaning simple vertical
inheritance doesn't explain their shared repertoire, and the
Fol4287–race3 SIX7/SIX10/SIX12 synteny result is the strongest evidence
in the dataset for accessory-region-level movement contributing to that
pattern.
