# Decision: IQ-TREE Version and Environment Choice

## Context

Two IQ-TREE installations were found available on the working machine:

- `iqtree` (version 3.1.3), installed inside the `six_project` conda environment.
- `iqtree2` (version 2.0.7), installed system-wide at `/usr/bin`, outside any conda environment.

A decision was required as to which installation to use for all phylogenetic tree-building steps in this project (core-genome tree and the six per-SIX-gene trees).

## Decision

The conda-environment installation, `iqtree` (version 3.1.3), was selected for all tree-building steps in this project.

## Rationale

**Reproducibility.** The project's environment is defined by the `six_project` conda environment, which is intended to be exportable and rebuildable independently of this specific machine. The system-wide `iqtree2` binary at `/usr/bin` falls outside that environment and would not be recreated by rebuilding `six_project` from an environment export. Depending on it would silently break reproducibility for anyone following the project's README, since the tool would appear installed on this machine but would be absent after a clean environment rebuild elsewhere.

**Version currency.** Version 3.1.3 is the more recent release relative to 2.0.7 and is the version actually integrated into the project's tracked environment, rather than a separate ad hoc system installation of unknown provenance.

## Consequence: command-line flag changes

IQ-TREE's flag names changed between the major version used in initial planning references (version 2.x conventions) and the version actually adopted (3.1.3). The following substitutions were confirmed directly against `iqtree --help` output before use, rather than assumed from prior documentation:

| Version 2.x convention | Version 3.1.3 confirmed equivalent | Purpose |
|---|---|---|
| `-bb NUM` | `-B NUM` (`--ufboot`) | Ultrafast bootstrap replicates |
| `-nt AUTO` | `-T AUTO` | Thread/core auto-detection |
| `-m MFP` | `-m MFP` (unchanged) | ModelFinder Plus: automatic best-fit substitution model selection followed by tree inference |
| `-seed NUM` | `--seed NUM` (unchanged in function) | Random seed for reproducible runs |

All tree-building commands in this project use the version 3.1.3 flag conventions listed above.

## Note on ultrafast bootstrap choice

Ultrafast bootstrap (`-B`, Minh et al. 2013) was used in preference to non-parametric bootstrap (`-b`) because it approximates branch support via resampling of estimated log-likelihoods (RELL) rather than fully re-optimizing the tree per replicate, making 1000 replicates computationally practical at the scale of this project (27 genomes; a core-genome supermatrix concatenated from multiple BUSCO loci). This is standard practice in current fungal phylogenomics and is the convention assumed by the project's reference materials.
