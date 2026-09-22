# Criteria for genuine SIX homolog identification

Percent-identity and query-coverage thresholds for accepting a candidate
locus (see docs/decisions/hsp_merging_criteria.md) as a genuine homolog of
a given SIX gene were set relative to documented sequence diversity within
each SIX gene family, rather than applying a single threshold uniformly
across all six genes.

Taylor et al. (2016, Molecular Plant Pathology) reported nucleotide
identity between Fusarium oxysporum f. sp. cepae and the f. sp. lycopersici
reference sequences ranging from 85% to 96% for SIX3, SIX5, SIX7, SIX10,
and SIX12, while sequences of these same genes were identical across all
isolates within a single forma specialis. This indicates that genuine,
distantly related orthologs of these five genes can legitimately fall well
below 100% identity when compared across formae speciales.

SIX1 shows greater documented variability than the above five genes.
Ghosal et al. (2024, Plant Pathology) reported that SIX1 protein length
varies both among formae speciales and among isolates within the same
forma specialis. In Fusarium oxysporum f. sp. cubense, SIX1 exists as
named allelic variants (SIX1a in race 1 versus SIX1b and SIX1c in tropical
race 4; PLOS ONE, 2018), and within F. oxysporum f. sp. lycopersici itself,
polymorphism in the related gene SIX3 is used to distinguish race 2 from
race 3 isolates (FEMS Microbiology Letters, 2009), indicating that
meaningful allelic variation can exist even within a single forma
specialis for this gene family.

On this basis, the following thresholds were applied to weighted percent
identity and query coverage (as calculated for each merged candidate
locus):

| Gene group | Identity floor | Confident identity | Coverage floor |
|---|---|---|---|
| SIX3, SIX5, SIX7, SIX10, SIX12 | 70% | 80% | 70% |
| SIX1 | 55% | 65% | Not enforced; noted but not used to reject |

A candidate locus with identity below the floor for its gene was classified
as not representing a genuine homolog. A candidate locus with identity at
or above the floor but below the confident threshold was retained but
flagged for manual inspection rather than automatically accepted, since
this range cannot be reliably distinguished from a spurious hit by
identity alone. For SIX3, SIX5, SIX7, SIX10, and SIX12, a candidate locus
meeting the identity threshold but falling below the coverage floor was
retained as a genuine homolog but flagged as a candidate for truncation,
pseudogenization, or assembly-boundary effects, to be resolved once the
underlying sequence is examined directly, rather than discarded. Coverage
was not used as a rejection criterion for SIX1, because the SIX1 reference
sequence used as the BLAST query is itself a partial, N-terminally
truncated coding sequence (see docs/decisions/protein_query_extraction.md),
so incomplete coverage against this reference does not necessarily
indicate an incomplete homolog in the target genome.
