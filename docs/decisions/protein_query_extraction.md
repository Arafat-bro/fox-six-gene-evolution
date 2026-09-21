# Protein and coding-sequence extraction for reference SIX genes

For each of the six reference accessions (SIX1: MK906592.1; SIX3: MK906598.1;
SIX5: MK906607.1; SIX7: GQ268954.1; SIX10: MK906667.1; SIX12: MW160867.1), the
protein query sequence used for homolog searching was taken directly from the
translation qualifier of the annotated coding-sequence (CDS) feature in the
corresponding NCBI GenBank record, rather than re-derived by translating the
raw nucleotide sequence. This qualifier reflects the original submitter's
coding-sequence interpretation, including the correct reading frame for
partial coding sequences and, for genes whose coding sequence spans multiple
exons, the correct removal of intervening intronic sequence prior to
translation.

The corresponding in-frame nucleotide coding sequence for each gene was
extracted using the CDS feature's annotated genomic location, which joins
exon coordinates in the correct order for genes with multi-part coding
sequences and applies reverse complementation where the feature is annotated
on the minus strand.

The SIX1 reference sequence is a partial, N-terminally truncated coding
sequence beginning mid-codon. As a result, percentage query coverage below
100% for SIX1 in subsequent homolog searches does not by itself indicate a
partial or degraded homolog in a target genome, since the reference query
itself does not represent a full-length protein.

The SIX5, SIX10, and SIX12 reference coding sequences contain introns.
Because translated nucleotide search operates on genomic (unspliced) target
sequence, a genuine homolog in a target assembly may be recovered as two or
more separate alignment segments corresponding to individual exons rather
than a single contiguous alignment. Such segments were merged prior to
determining overall percentage identity and query coverage for a candidate
homolog.
