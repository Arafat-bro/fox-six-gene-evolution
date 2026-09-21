## Genome Selection

Twenty-seven *Fusarium oxysporum* genome assemblies were selected from
NCBI GenBank to serve as the comparative dataset for this analysis.
Candidate assemblies were retrieved from the full set of publicly
available *F. oxysporum* genome records using the NCBI Datasets
command-line tool, restricted to the GenBank collection rather than
RefSeq, as GenBank contains substantially greater representation of the
host and forma specialis diversity required for this study.

Assemblies were prioritized for selection based on assembly contiguity
(contig N50) and completeness of associated metadata, with Chromosome-
level or Complete Genome assemblies preferred where available. Where
multiple candidate genomes existed for the same host, forma specialis,
or race, redundant or near-identical isolates were excluded in favor of
higher-contiguity representatives, to avoid overrepresenting a single
population at the expense of broader diversity. Conversely, isolates
representing genuinely distinct lineages, races, or geographic origins
within the same forma specialis were retained even where this meant
including more than one genome per host, to support within-group as
well as between-group comparison.

Two strains were required by design and included regardless of the
general selection criteria: *F. oxysporum* f. sp. *lycopersici* strain
Fol4287, the reference strain in which the SIX7-SIX10-SIX12
pathogenicity chromosome region was originally characterized, and Fo47,
a nonpathogenic isolate included as a negative control. For Fol4287, an
improved assembly of the reference strain (GCA_003315725.1) was used in
place of the original 2010 assembly (GCA_000149955.2), owing to
substantially higher contig N50 (1.34 Mb versus 95 kb) and correction of
sequencing errors present in the original low-coverage assembly.

Forma specialis assignment for each genome was determined primarily from
the assembly's own submitted NCBI taxonomy record. Where this
designation was not present in a given assembly's record, forma
specialis was instead assigned based on established host-association
literature; these cases are noted in the genome metadata table. One
notable case, an isolate collected from cowpea (*Vigna unguiculata*,
GCA_054643655.1), is submitter-designated *F. oxysporum* f. sp.
*phaseoli* (the bean-associated lineage) rather than f. sp.
*tracheiphilum*, the lineage typically associated with cowpea; this
designation was retained as submitted and is discussed further in the
interpretation of host-association patterns.

One candidate genome isolated from a human clinical case was excluded,
as this analysis is restricted to plant-associated isolates consistent
with its focus on host-associated pathogenicity.

The final dataset comprises 27 genomes spanning 15 distinct hosts and
associated formae speciales, including paired pathogenic and
nonpathogenic isolates from two independent host systems (tomato:
Fol4287/Fo47; koa: *F. oxysporum* f. sp. *koae*/a nonpathogenic isolate
from the same study). Full genome-level metadata, including strain,
host, geographic origin, assembly statistics, and sequencing technology,
is provided in the accompanying genome metadata table
(`results/genome_metadata.tsv`).
