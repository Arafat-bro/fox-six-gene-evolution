# Criteria for merging BLAST high-scoring segment pairs into candidate homolog loci

Because three of the six reference SIX genes (SIX5, SIX10, SIX12) contain
introns, a single genomic copy of these genes in a target assembly is
frequently returned by translated nucleotide search (tblastn) as multiple
separate high-scoring segment pairs (HSPs), one per exon, rather than a
single contiguous alignment.

Two HSPs found on the same scaffold and strand were merged into a single
candidate homolog locus (interpreted as one gene copy) if both of the
following conditions were met:

1. The genomic distance between the two HSPs did not exceed 3000 base
   pairs. This threshold was chosen to comfortably exceed the expected
   size of introns in these effector genes while remaining well below the
   genomic distance expected between independent gene copies.
2. The query (reference protein) coordinate ranges of the two HSPs
   overlapped by no more than 15 amino acids. Genuine exon fragments of a
   single gene are expected to tile across the length of the reference
   protein with minimal or no overlap in the region of the query each HSP
   aligns to; extensive overlap instead indicates that two HSPs align to a
   similar portion of the reference protein and are therefore more
   consistent with two independent, closely spaced gene copies than with
   exon fragments of one gene.

HSP pairs that were genomically close but exceeded the query-overlap
tolerance were not merged automatically. These cases were instead recorded
separately for manual inspection, consistent with the requirement to
distinguish automated results from conclusions based on manual review.

For each resulting candidate locus, overall query coverage was calculated
as the union of the query positions covered by all merged HSPs (avoiding
double-counting of any overlapping residues) divided by the total length of
the reference protein. Overall percent identity was calculated as the
length-weighted mean identity across all merged HSPs.
