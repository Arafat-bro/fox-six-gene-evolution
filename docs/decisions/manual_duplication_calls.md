## Manual resolution of ambiguous scaffold-duplication calls

Two (gene, genome) pairs were flagged by the automated deduplication script
(`05_dedupe_scaffold_duplicates.py`) as `assembly_duplicate_ambiguous_which_to_keep`,
because both candidate loci in each pair sat on non-chromosome-level scaffolds,
so the script's automated tie-breaker (prefer the chromosome-assigned accession)
could not be applied. Both cases were investigated manually by direct coordinate
inspection and sequence extraction, rather than relying on BLAST summary statistics
alone.

### Case 1: SIX3 / GCA_003615085.1 (copy1 vs copy2, match_type = identical)

Copy 1 was located at MRCU01000016.1:162130-162615 (minus strand); copy 2 at
MRCU01000019.1:216560-217045 (minus strand). Both hits reported identical BLAST
statistics (percent identity, alignment length, coverage, e-value, and bitscore)
to the decimal place, which is not expected of two independently evolving loci
and was treated as grounds for direct verification.

The underlying assembly (FoC_Fus2, PacBio + Illumina, 34 contigs, N50 4.14 Mb)
is chromosome-level overall, but both candidate loci fall on small unplaced
scaffolds (MRCU0100####) rather than on assigned chromosome pseudomolecules --
consistent with the known failure mode in which a small contig is not fully
collapsed during assembly and its sequence is represented twice under separate
accessions.

Both loci were extracted directly from the genome FASTA with `seqkit subseq`,
reverse-complemented to a common orientation, and compared directly. The
sequence bodies were confirmed identical in full (486 bp); only the FASTA
headers, which necessarily differ by scaffold accession and coordinates, showed
any difference.

**Conclusion:** the two loci represent a single physical copy of SIX3 duplicated
across two scaffold accessions by an assembly artifact, not two independent gene
copies. One copy is excluded from the deduped homolog table.

**Which copy was retained:** because the two sequences are identical, retaining
either copy has no effect on any downstream sequence-based analysis (alignment,
tree-building, dN/dS). The automated script's existing tie-breaker (prefer the
chromosome-level accession) does not apply, as neither scaffold is chromosome-
assigned. In the absence of an existing precedent for this situation, the lower-
numbered scaffold accession was retained (MRCU01000016.1; copy2 on MRCU01000019.1
excluded) as a simple, consistent, and explicitly arbitrary convention, to be
applied to any future ties among non-chromosome-level scaffolds for consistency.

### Case 2: SIX7 / GCA_013423245.1 (copy2 vs copy4, match_type = reverse_complement)

Both loci lie on the same scaffold, CM023981.1, a chromosome-level accession,
but approximately 3.7 Mb apart (positions ~304,130-304,771 and
~4,043,030-4,043,671), on opposite strands.

**Conclusion:** this separation cannot be explained by an uncollapsed-contig
duplication artifact, which produces overlapping or immediately adjacent
coordinates, not multi-megabase separation on a single finished chromosome.
The two loci are therefore treated as genuinely independent gene copies rather
than a duplicate scaffold representation of one locus.

**Decision:** both copies are retained in the deduped homolog table. This genome
carries a four-copy SIX7 expansion overall and is flagged as a candidate case of
lineage-specific gene duplication, to be carried forward into the Part 9
integration table alongside the phylogenetic and presence/absence evidence.

### Case 3: SIX1 / GCA_016164145.2 (copy2 vs copy3, borderline bucket)

Copy 2 was located at WILW02000022.1:71194-72018 (plus strand); copy 3 at
CM040507.1:2449316-2450140 (minus strand). Both loci reported identical
BLAST statistics (100% query coverage, 64.06% identity, matching e-value and
bitscore pattern differing only by strand), and both fell into the
borderline identity bucket rather than the genuine bucket.

This case differs in kind from Cases 1 and 2 above: it was not flagged as
ambiguous by an automated script, because neither automated
scaffold-duplication check ever inspected it. `05_dedupe_scaffold_duplicates.py`
compares loci only within the genuine bucket; `06_apply_manual_dedup_calls.py`
reads only from the genuine-bucket dedup output. A duplicate pair confined
entirely to the borderline bucket therefore fell through a genuine scope
gap in the automated pipeline, rather than being evaluated and left
ambiguous.

Both loci were extracted directly with seqkit subseq, oriented to a common
strand (copy 3 reverse-complemented), and compared directly. The sequence
bodies were confirmed identical in full (~825 bp); only the FASTA headers
differed.

**Conclusion:** the two loci represent a single physical copy of SIX1
duplicated across two scaffold accessions by an assembly artifact (the same
uncollapsed-haplotig failure mode as Case 1), not two independent gene
copies. Copy 2 is excluded; copy 3 is retained, consistent with the
lower-numbered-accession convention is not applicable here (WILW vs CM
prefixes are not comparable under that rule), so instead the
chromosome-style-accession preference used by the automated dedup script's
own tie-breaker was applied: CM040507.1 is a chromosome-level accession
(CM prefix) while WILW02000022.1 is a WGS contig-style accession, so copy 3
(CM040507.1) was retained as the chromosome-assigned representative.

Note on provenance: this duplicate pair was first suggested by output from
a separate AI assistant session with no access to this project's actual
files or verified history. Its claimed file path for the relevant script
was independently checked and found incorrect, and one of its stated
supporting details (a specific mismatch count over a stated window size)
did not match this session's own direct verification. The core claim of a
duplicate pair at these two loci was independently confirmed against the
project's own primary sequence data before being acted upon; no detail from
that external session was taken on trust.
