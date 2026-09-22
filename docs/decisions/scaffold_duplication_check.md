# Detection of assembly-duplicated loci

Two candidate SIX1 loci identified in genome GCA_016164145.2
(SIX1_GCA_016164145.2_copy2, on scaffold WILW02000022.1; and
SIX1_GCA_016164145.2_copy3, on scaffold CM040507.1) showed identical
percent identity, e-value, and bitscore against the SIX1 reference
sequence despite residing on different scaffolds and opposite strands.
Direct extraction and comparison of the underlying genomic sequence at
each locus showed the two sequences to be exact reverse complements of
one another (825 bp each). This indicates that a single genomic region
was represented twice in this draft assembly, rather than representing
two independent gene copies.

For copy-number purposes, this genome was treated as carrying a single
SIX1 copy at this locus. The scaffold WILW02000022.1 hit was excluded as
a duplicate assembly representation; the CM040507.1 hit was retained, as
CM040507.1 follows a chromosome-level naming convention in this assembly
whereas WILW02000022.1 follows a WGS contig naming convention, indicating
CM040507.1 is the more reliable representation of this genomic region.

All 50 candidate loci surviving identity/coverage filtering were checked
programmatically for this same pattern (same gene, same genome, multiple
loci) prior to finalizing the presence/absence and copy-number matrix; see
scripts/03_six_homolog_search/05_dedupe_scaffold_duplicates.py.
