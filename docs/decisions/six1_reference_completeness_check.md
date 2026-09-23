## Confirmation: SIX1 reference protein correctly handles its incomplete GenBank annotation

The SIX1 reference GenBank record (MK906592.1) annotates its CDS as
`<1..845` with `codon_start=3`. The `<1` indicates the deposited sequence is
incomplete at the 5' end (missing an unknown number of residues from the
true N-terminus of the protein), and `codon_start=3` indicates that, within
this incomplete record, translation must begin at position 3 rather than
position 1 to remain in the correct reading frame.

This was checked to confirm that `SIX1_protein.fasta`
(data/references/processed/), used throughout this project as the tblastn
query and as the reference-length anchor for truncation classification, was
built with this frame offset correctly applied, rather than by naively
translating the raw nucleotide record from position 1.

Verification: the first several dozen residues of SIX1_protein.fasta were
compared directly against the `/translation=` qualifier in the GenBank
record itself (which NCBI computes with codon_start already correctly
applied). The two matched exactly. This confirms the reference protein
query was extracted correctly with respect to reading frame; no correction
is required.

Remaining limitation, not a pipeline error: because the GenBank record
itself is annotated as N-terminally incomplete, the reference protein used
throughout this project is itself somewhat shorter than the true full-length
SIX1 protein. Any locus classified as "truncated" for SIX1 (protein length
below 90% of the reference protein's length, per
04_homolog_validation/01_extract_translate_classify.py) is being compared
against a reference that is already known to be incomplete at its N-terminus.
This is a property of the deposited reference record and cannot be corrected
without an alternative, complete SIX1 reference sequence; it is noted here
so it can be stated explicitly as a limitation in the project's methods or
limitations section, rather than surfacing as an unexplained inconsistency.
