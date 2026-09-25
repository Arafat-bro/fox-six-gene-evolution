# Foxy_SIX_Evolution — SIX Effector Gene Evolution in Fusarium oxysporum

## Project structure and script order
1. Genome selection + metadata: results/genome_metadata.tsv (27 genome assemblies, manually curated from NCBI)
2. Core-genome marker identification: scripts/05_core_genome/miniprot_single_gene/01_make_queries.py
   - BUSCO (fungi_odb10) on 6 representative genomes → 3,912 shared single-copy candidate genes
3. Core-genome ortholog mapping: scripts/05_core_genome/miniprot_single_gene/02_run_miniprot.sh
   - Miniprot v0.18 protein-to-genome alignment against all 27 genomes
4. Core-genome ortholog filtering: scripts/05_core_genome/miniprot_single_gene/03_extract.py
   - MINCOV=0.7, MINFRAC=0.9 → 3,911 of 3,912 genes retained
5. Core-genome alignment/tree: scripts/05_core_genome/miniprot_single_gene/04_align_tree.sh
   - MAFFT --auto, trimAl -automated1, concatenation, IQ-TREE 3 (LG+I+G, 1000 UFBoot)
   - Output: results/core/core_tree_fast.treefile
6. SIX gene homolog search: tblastn of 6 reference SIX genes against all 27 genomes
   - Thresholds: e-value ≤1e-10, query coverage ≥70%
   - Output: results/six_presence_absence_matrix.tsv, results/six_pseudogene_annotation_matrix.tsv

## Environment
conda environment `six_project` — see environment.yml for pinned versions
(Python 3.13.15, Miniprot 0.18, BUSCO 6.0.0, MAFFT 7.526, trimAl 1.5.1, IQ-TREE 3.1.3)

## Manual steps
- Genome selection criteria and forma specialis/host/geography metadata were curated manually
  from NCBI BioSample records (results/genome_metadata.tsv)
- Borderline BLAST hits were flagged for manual review rather than automatically included/excluded

## Not completed (see Limitations in report)
- Per-SIX-gene phylogenies and tree-vs-core-tree congruence comparison
- dN/dS (HyPhy FEL) analysis
- SIX7–SIX10–SIX12 synteny/linkage analysis (clinker)
