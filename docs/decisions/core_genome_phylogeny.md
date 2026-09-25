# Core-genome phylogeny: methods and decisions

## Problem
BUSCO (hypocreales_odb12.2, `--miniprot`) took about 1 hour per genome and
crashed twice when run in parallel on an 8 GB RAM / 8 core machine. Full
BUSCO on all 27 genomes was not feasible in the time available. Only 6
genomes had finished BUSCO runs:
GCA_003315725.1, GCA_003615085.1, GCA_013085055.1, GCA_049306905.1,
GCA_053754965.2, GCA_055944325.1.

## Decision
Reuse the single-copy BUSCO protein set from the 6 finished genomes as
query proteins, and map them directly to all 27 genome assemblies with
`miniprot`, bypassing BUSCO's full pipeline for the remaining 21 genomes.
This was chosen over:
- Installing `compleasm` (dependency conflict with six_project's
  Python 3.13 environment; also may not support odb12).
- A Mash/fastANI distance tree (rejected as the primary tree since it
  gives no gene alignment, no substitution model, and no branch support;
  kept only as a fallback that was not ultimately needed).

## Pipeline
1. `01_make_queries.py` — read the single-copy BUSCO protein sets from
   the 6 finished genomes, take the genes shared by all 6, and write one
   reference protein per gene (from GCA_003615085.1) to
   `results/core/queries.faa`.
2. `02_run_miniprot.sh` — run `miniprot -t 8 --gff -N 0` for each of the
   27 genome FASTAs in `data/genomes/raw/flat/` against `queries.faa`,
   one genome at a time (parallel runs previously crashed the machine).
   Writes to a `.gff.tmp` then renames to `.gff` so an interrupted run is
   retried rather than treated as complete.
3. `03_extract.py` — parse each GFF, keep the best-scoring hit per query
   gene, require query coverage >= 0.7, extract and translate the CDS,
   and keep only genes present in at least 90% of the 27 genomes.
   Produced 3,911 shared genes.
4. `04_align_tree.sh` — align each gene with MAFFT (`--auto`), trim with
   trimAl (`-automated1`), concatenate per genome
   (`results/core/concat.faa`, ~2.3M columns), and build the tree.

## Deviations from the original plan (and why)
- The concatenated alignment (~2.3M columns) was subsampled to ~250k
  columns before tree inference, and `iqtree3` was run with `--fast` and
  `-alrt 1000` (SH-aLRT) rather than full ModelFinder + 1000 UFBoot, to
  fit the time constraint. Model used: LG+I+G.
- Final tree: `results/core/core_tree_final.treefile`
  (LogLikelihood = -1,211,567.3, tree length = 0.252).
- This is a single-model, subsampled, SH-aLRT-supported tree rather than
  a fully partitioned, full-length, UFBoot-supported tree. This is a
  methodological limitation to state explicitly in the report.

## Known caveats
- Query proteins come from one genome (GCA_003615085.1), so only genes
  that genome carries in single copy could be used as queries.
- Coverage/identity filters catch the worst frameshift or fragment
  artifacts but were not manually reviewed hit-by-hit; treat as
  automated results, not manually curated ones.
