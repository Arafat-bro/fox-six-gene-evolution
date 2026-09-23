# Decision: BUSCO Tool and Lineage Dataset Choice

## Context

The assignment specifies building a core-genome phylogeny from conserved single-copy orthologs but does not mandate a specific tool for identifying them. BUSCO (`hypocreales_odb12.2`) was selected for this step.

## Decision

BUSCO was chosen over OrthoFinder because it works directly on raw, unannotated genome assemblies and draws on a pre-curated, taxonomy-specific single-copy marker set, avoiding the need to first generate gene annotations for all 27 genomes. This is a standard, well-established approach for this kind of within-species phylogenomic backbone tree.

Within BUSCO's lineage datasets, `hypocreales_odb12.2` (4077 markers) was chosen over the broader `fungi_odb12.2` (1019 markers) because all 27 genomes belong to a single species within Hypocreales; a narrower, taxonomically appropriate lineage dataset yields substantially more usable single-copy markers with no loss of applicability, directly improving the resolution of the resulting supermatrix and tree.

Note: this project's reference materials name `fungi_odb10`; that dataset is no longer offered under the installed BUSCO version (6.0.0), which uses OrthoDB v12 (`_odb12.2`) datasets. `hypocreales_odb12.2` is the actual dataset used and reported here.
