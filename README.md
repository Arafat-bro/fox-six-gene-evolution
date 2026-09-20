# Fusarium oxysporum SIX Effector Gene Evolution

Comparative genomics analysis of SIX effector gene (SIX1, SIX3, SIX5, SIX7, SIX10, SIX12) distribution and evolutionary history across *Fusarium oxysporum* genomes.

## Environment setup

```bash
mamba env create -f environment.yml
conda activate six_project
```

## Folder structure

- `data/references/` — reference SIX gene sequences (raw NCBI downloads + processed)
- `data/genomes/` — target genome assemblies and their metadata
- `data/annotations/` — genome annotation files, where available
- `results/` — pipeline outputs (homolog tables, matrices, trees, synteny, selection stats)
- `scripts/` — all analysis scripts, numbered by pipeline stage
- `docs/decisions/` — manual judgment calls and their justification (genome selection criteria, threshold choices, verification findings)
- `docs/methods/` — detailed methods notes per pipeline stage
- `figures/` — final figures for the report
- `logs/` — captured output from long-running commands

## Reproducing this analysis

| Script | Produces | Notes |
|---|---|---|
| `scripts/01_download_references/download_reference_genes.sh` | `data/references/raw/*.fasta`, `*.gb` | Downloads the six NCBI reference accessions |
| `scripts/01_download_references/verify_reference_genes.sh` | terminal output only (manually reviewed) | Findings recorded in `docs/decisions/reference_gene_verification.md` |

*(this table grows as each pipeline stage is completed)*

## Manual judgment calls

See `docs/decisions/` for every point where a human decision was made (genome selection criteria, BLAST threshold choices, borderline hit calls, synteny interpretation) rather than a fully automated pipeline step.
