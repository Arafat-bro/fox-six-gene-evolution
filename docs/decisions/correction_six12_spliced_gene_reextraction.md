## Correction: SIX12 loci affected by the same intron-spanning frameshift issue as SIX5/SIX10, via a different mechanism

Following the SIX5/SIX10 correction (see
correction_six10_six5_spliced_gene_reextraction.md), SIX3 and SIX12's
reference GenBank records were also checked for splicing, which had not
previously been done. SIX12 (MW160867.1) is confirmed spliced
(`CDS join(1..22,71..432)`), a two-exon gene with a short (~48 bp) intron.

Four of five candidate SIX12 loci showed an identical, gene-specific
signature (protein_length=143 against a reference length of 127, with
exactly one internal stop codon) across independent, unrelated genomes --
the same "too consistent to be coincidence" pattern that originally
flagged the SIX5/SIX10 bug. However, unlike the SIX5/SIX10 cases, all four
affected SIX12 loci had n_hsps_merged=1, meaning they were not caught by
the diagnostic used previously (checking for merged multi-HSP loci). This
is because SIX12's intron is short enough to be absorbed as an internal
gap within a single reported BLAST HSP, rather than splitting the
alignment into two separate HSPs as occurred for SIX5 and SIX10's larger
introns. The one unaffected SIX12 locus (GCA_003615085.1) evidently did
not have its alignment span the intron.

Correction: the four affected loci were re-extracted using the same
exonerate (protein2genome) approach used for SIX5/SIX10, correctly
excising the intron before translation. All four resolved to complete,
full-length translations with zero internal stop codons.

Methodological note for any future gene: neither "check n_hsps_merged>1"
nor "check for a join() CDS" alone is a sufficient screen. A gene can be
genuinely spliced (per its GenBank record) yet still produce single-HSP
loci for a subset of genomes, if its intron is short enough for BLAST to
bridge internally. Reference splicing status must be checked directly
against the GenBank record for every gene, and classification output
should always be screened for anomalously identical stats across unrelated
genomes as a general practice, regardless of n_hsps_merged.
