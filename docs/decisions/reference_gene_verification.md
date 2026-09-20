# Reference SIX Gene Verification Notes

Verification performed on the six NCBI reference accessions before use as BLAST query sequences, to confirm CDS boundaries, completeness, and intron structure.

| Gene  | Accession   | CDS span              | Completeness | Introns | codon_start | Notes |
|-------|-------------|------------------------|---------------|---------|-------------|-------|
| SIX1  | MK906592.1  | <1..845                | Partial (5' truncated) | No | 3 | Query itself lacks true start; expect <100% coverage on genuinely full-length homologs. Translate from position 3, not 1. |
| SIX3  | MK906598.1  | 80..571                | Complete       | No | 1 | Clean single-exon CDS. |
| SIX5  | MK906607.1  | join(57..100,166..275,328..395,445..582) | Complete | Yes (4 exons) | 1 | Must merge multi-HSP BLAST hits before judging coverage/completeness. |
| SIX7  | GQ268954.1  | 1..663                 | Complete       | No | 1 | Clean single-exon CDS. |
| SIX10 | MK906667.1  | join(21..100,171..540) | Complete       | Yes (2 exons) | 1 | Same multi-HSP merging caveat as SIX5. |
| SIX12 | MW160867.1  | join(1..22,71..432)    | Complete       | Yes (2 exons) | 1 | Same multi-HSP merging caveat as SIX5/SIX10. |

## Implications for downstream pipeline

1. **Protein query extraction for BLAST searches:** for spliced genes (SIX5, SIX10, SIX12), the protein query used for `tblastn` must come from the *spliced* coding sequence (translate the joined exons together), not a straight substring of the genomic FASTA — a straight substring would include intron sequence and produce a garbage/frameshifted translation. Biopython's `SeqFeature.location.extract()` on the CDS feature handles `join()` automatically and is the safe way to do this.
2. **BLAST hit interpretation:** for SIX5/SIX10/SIX12, a single genuine homolog in a target genome may appear as 2+ separate HSPs in the `tblastn` output if that genome's copy also retains introns in the same positions. Before applying the coverage/identity thresholds, HSPs from the same (query, subject) pair that are colinear and closely spaced must be merged, and total coverage computed across the merged span — not per individual HSP row.
3. **SIX1 coverage caveat:** because the reference query itself is partial (missing the 5' end), sub-100% query coverage on SIX1 hits should not be interpreted as evidence of truncation/pseudogenization in the target genome without checking the extracted sequence by eye first.
