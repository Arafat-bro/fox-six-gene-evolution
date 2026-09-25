## Reference SIX Gene Sequences

Reference sequences for the six SIX effector genes examined in this
study — SIX1, SIX3, SIX5, SIX7, SIX10, and SIX12 — were retrieved from
NCBI GenBank using the accessions MK906592.1, MK906598.1, MK906607.1,
GQ268954.1, MK906667.1, and MW160867.1 respectively. For each gene, both
the nucleotide FASTA sequence and the full GenBank flatfile record were
downloaded via NCBI's Entrez utilities (`efetch`), rather than the
FASTA sequence alone, because the flatfile's CDS feature annotation is
required to establish each gene's coding boundaries, reading frame, and
intron/exon structure before the sequence can be used as a reliable
BLAST query or reference length for downstream classification.

Each reference record's CDS annotation was inspected individually prior
to use, rather than assumed uniform across all six genes, since GenBank
records can vary in completeness and in whether they represent spliced
mRNA or unspliced genomic sequence. This inspection identified three
spliced reference genes, one reference gene with an incomplete record,
and two straightforward unspliced, complete records:

**Spliced (multi-exon) reference genes.** SIX5 (MK906607.1, four exons),
SIX10 (MK906667.1, two exons), and SIX12 (MW160867.1, two exons, with a
short ~48 bp intron) are each annotated with multi-segment `join()` CDS
coordinates, confirming all three are intron-containing genes represented
as spliced mRNA. This distinction proved essential in practice: candidate
genomic loci for these genes were in some cases detected by BLAST as
multiple separate high-scoring segment pairs (one per exon) and in other
cases, where the intron was short, as a single alignment spanning the
intron internally. In both cases the intron sequence must be identified
and excluded before translation. An initial extraction method that
treated each locus as a single contiguous genomic span was found to
introduce frameshift errors at intron boundaries, producing translations
that falsely resembled pseudogenization, for both multi-segment loci
(SIX5, SIX10) and, less obviously, for single-segment loci whose
alignment happened to bridge a short intron (SIX12). All affected loci
across all three genes were identified — the latter case specifically by
the presence of identical, gene-specific translation statistics across
otherwise-unrelated genomes — and corrected using exonerate
(protein2genome model) to identify true exon boundaries prior to
translation. Full detail is provided in
`docs/decisions/correction_six10_six5_spliced_gene_reextraction.md` and
`docs/decisions/correction_six12_spliced_gene_reextraction.md`.

**Unspliced reference gene with an incomplete record.** SIX1
(MK906592.1) is annotated with a partial CDS (`<1..845`,
`codon_start=3`), indicating that the deposited record is missing an
unknown number of residues from the true N-terminus of the protein,
and that translation of this particular record must begin at its third
nucleotide position to remain in the correct reading frame. The
reference protein sequence used throughout this study was confirmed to
match NCBI's own precomputed translation of this record exactly,
confirming the frame offset had already been correctly applied at
extraction. This incompleteness is a property of the deposited record
itself, not a correctable error, and is treated as a limitation
(`docs/decisions/six1_reference_completeness_check.md`): loci classified
as truncated relative to this particular reference should be understood
as truncated relative to an already N-terminally incomplete sequence, not
the full-length native protein.

**Unspliced, complete reference genes.** SIX3 (MK906598.1, `CDS 80..571`)
and SIX7 (GQ268954.1, `CDS 1..663`) are each annotated as single, complete
coding sequences with no splicing and no frame-offset complication.

For each reference gene, the corresponding protein sequence was
extracted for use as the tblastn query, using NCBI's own precomputed
translation from the GenBank record (accounting for `codon_start` and
completeness flags as above) rather than translating the raw nucleotide
sequence independently, to avoid introducing a frame or boundary error
at the extraction step itself.