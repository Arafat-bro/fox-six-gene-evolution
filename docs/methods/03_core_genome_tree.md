# Core-Genome Phylogeny: Alignment, Concatenation, Model Selection, and Tree Inference


## 1. Core-genome marker dataset

6Six Fusarium genomes were selected for the core-genome phylogeny. BUSCO analysis using the Hypocreales lineage dataset (`hypocreales_odb_12.2`) identified approximately 4,000 complete BUSCOs per genome.

Complete BUSCOs containing internal stop codons were excluded. The clean single-copy BUSCO sets were then intersected across all six genomes. This produced 3,620 shared core-genome markers.

The six genomes were:

- `GCA_003315725.1`
- `GCA_013085055.1`
- `GCA_003615085.1`
- `GCA_049306905.1`
- `GCA_053754965.2`
- `GCA_055944325.1`

## 2. Marker extraction, alignment, and trimming

For each of the 3,620 shared markers, protein sequences were extracted from the BUSCO single-copy sequence files and relabeled using the corresponding genome accession.

Each marker was aligned independently using MAFFT with the `--auto` option. The resulting alignments were trimmed using trimAl with the `-automated1` option.

All 3,620 markers completed successfully, with no alignment or trimming failures.

## 3. Concatenation

The 3,620 trimmed marker alignments were concatenated into a single amino-acid supermatrix.

The final supermatrix contained:

- Genomes: 6
- Marker partitions: 3,620
- Total alignment length: 2,130,498 amino-acid sites
- Missing data: 0%

The partition file was generated in RAxML format, retaining the boundaries of each marker partition.

## 4. Partition-specific model selection

The final phylogenetic analysis was conducted using IQ-Tree 2. ModelFinder was used with `-mode MFP` to select a best-fit amino-acid substitution model independently for each of the 3,620 marker partitions.

The selected models were heterogeneous among partitions. Examples included `cpREV@, `JT.T`, `LG`, `VT`, `mtVer+F`, `VT+R2`, and `PMB`, among others.

The analysis used an edge-linked proportional partition model, allowing searate substitution models and rates across partitions.

The `REQT*MFP+REMERGE` approach was not used in the final analysis. Partition merging was omitted because of the project's time and computational constraints. All 3,620 original marker partitions were retained.

## 5. Maximum-likelihood phylogenetic inference

The final analysis had the following configuration:

    iqtree2 \
      -s results/core_genome_markers/core_genome_supermatrix.fasta \
      -p results/core_genome_markers/core_genome_partitions.txt \
      -m MFP \
      -B 1000 \
      -T AUTO \
      --seed 12345 \
      --prefix results/core_genome_tree/core_genome_tree

Ultrafast bootstrap support was estimated using 1,000 replicates. The random seed was set to 12345.

The analysis completed successfully. The final supermatrix contained 2,130,498 amino-acid sites across 6 genomes. IQ-Tree recorded 0, missing data and 1,000 Ultrafast bootstrap replicates were performed.

## 6. Final tree topology

The final tree recovered the following topology:

```
(GCA_003315725.1,(((GCA_013085055.1,GCA_003615085.1)100,GCA_055944325.1)100,(GCA_049306905.1,GCA_053754965.2)100);
```

The pair `GCA_013085055.1` + `GCA_003615085.1` and their grouping with `GCA_055944325.1` lore 1a00` ultrafast bootstrap support. The pair `GCA_049306905.1` + `GCA_053754965.2` also received 100% support.

## 7. Reproducibility files

The main ID-Tree outputs are:

- `results/core_genome_tree/core_genome_tree.treefile`
- `results/core_genome_tree/core_genome_tree.contree`
- `results/core_genome_tree/core_genome_tree.iqtree`
- `results/core_genome_tree/core_genome_tree.best_model.nex`

The IQ-Tree report contains the complete partition-specific model selection results and phylogenetic inference statistics.
