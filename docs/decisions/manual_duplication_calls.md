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
