# Genome Selection Criteria and Manual Judgment Calls

## Host scope
This analysis is restricted to *Fusarium oxysporum* isolates recovered from
plant hosts. One candidate genome (GCA_009746015.2) was isolated from a
human clinical case and was excluded, as clinical isolates fall outside
the scope of a plant-pathogenicity/host-association analysis.

## Required strain inclusions
- **Fol4287** (GCA_003315725.1): an improved assembly of the original
  reference strain in which the SIX7-SIX10-SIX12 pathogenicity chromosome
  region was first characterized, selected over the original 2010 assembly
  (GCA_000149955.2) for substantially higher contig N50 (1.34 Mb vs. 95 kb)
  and correction of sequencing errors present in the original low-coverage
  assembly.
- **Fo47** (GCA_013085055.1): included as the required nonpathogenic
  comparison strain. As a nonpathogenic isolate, it has no associated host
  in NCBI's BioSample record and was therefore added manually rather than
  emerging from the host-availability filter applied to pathogenic isolates.

## Banana (Fusarium oxysporum f. sp. cubense) isolate selection

Five candidate banana-associated genomes passed initial quality and host
availability filtering:

- GCA_011316005.3 (TC1-1): Foc Race 1 (Thangavelu et al. 2021, *Comparative
  Whole-Genome Sequence Analyses of Fusarium Wilt Pathogen (Foc R1, STR4
  and TR4) Infecting Cavendish (AAA) Bananas in India*, J Fungi)
- GCA_014282265.3 (BC2-4): Foc TR4 (Thangavelu et al. 2021, ibid.;
  Thangavelu et al. 2022, *Development of PCR-Based Race-Specific Markers
  for Differentiation of Indian Fusarium oxysporum f. sp. cubense*, J Fungi)
- GCA_031834405.1 (II5): Foc TR4 (Berg et al. 2012)
- GCA_053754965.2 (Tw-TR4-1): Foc TR4
- GCA_053754955.2 (Tw-R1-1): Foc Race 1

TC1-1 and BC2-4 were excluded from the final dataset. Both have markedly
lower contig N50 (95 kb and 18 kb, respectively) than the other available
genomes of the same race, and each race is already represented in the
dataset by a substantially more contiguous assembly (Tw-R1-1 for Race 1;
II5 and Tw-TR4-1 for TR4). The final dataset therefore retains three
banana genomes spanning two races (Race 1, TR4), each with contig N50
above 4.4 Mb.

## Duplicate and redundant assembly removal

- GCA_014154955.1 (Brassica oleracea isolate Fo5176, Chromosome-level,
  N50 3.38 Mb) was removed as a duplicate of GCA_030345115.2, a Complete
  Genome assembly of the same isolate (Fo5176) with higher contiguity
  (N50 4.37 Mb). Only the higher-quality assembly was retained.
- GCA_027920445.1 (Fusarium odoratissimum/oxysporum f. sp. cubense,
  strain hn51, NCBI BioProject PRJNA918855) was excluded. This isolate
  is confirmed Tropical Race 4, a lineage already represented in the
  dataset by two independently sequenced isolates (Tw-TR4-1, II5); its
  inclusion would add a third near-identical TR4 representative without
  contributing additional race-level diversity.

## Acacia koa isolate pair

Two Acacia koa-associated genomes were identified as a deliberately paired
pathogenic/nonpathogenic isolate set from the same study (NCBI BioSample
accessions SAMN14669697 and SAMN14669696), rather than duplicate assemblies
of a single isolate:

- GCA_014857105.1 (isolate "Fo koae 44"): Fusarium oxysporum f. sp. koae,
  the pathogenic isolate. Chromosome-level, 310 contigs, N50 4.22 Mb.
- GCA_014857085.1 (isolate "170"): Fusarium oxysporum, no forma specialis
  assigned; the nonpathogenic comparison isolate from the same study.
  Chromosome-level, 273 contigs, N50 587 kb.

Both genomes were retained. This pair provides an independent
pathogenic/nonpathogenic comparison outside of the Fol4287/Fo47 pair,
supporting cross-validation of any observed association between SIX gene
presence and pathogenicity. The lower contiguity of GCA_014857085.1 is
noted as a caveat for gene-absence calls specific to that genome.

## Acacia koa isolate pair

Two Acacia koa-associated genomes were identified as a deliberately paired
pathogenic/nonpathogenic isolate set from the same study (NCBI BioSample
accessions SAMN14669697 and SAMN14669696), rather than duplicate assemblies
of a single isolate:

- GCA_014857105.1 (isolate "Fo koae 44"): Fusarium oxysporum f. sp. koae,
  the pathogenic isolate. Chromosome-level, 310 contigs, N50 4.22 Mb.
- GCA_014857085.1 (isolate "170"): Fusarium oxysporum, no forma specialis
  assigned; the nonpathogenic comparison isolate from the same study.
  Chromosome-level, 273 contigs, N50 587 kb.

Both genomes were retained. This pair provides an independent
pathogenic/nonpathogenic comparison outside of the Fol4287/Fo47 pair,
supporting cross-validation of any observed association between SIX gene
presence and pathogenicity. The lower contiguity of GCA_014857085.1 is
noted as a caveat for gene-absence calls specific to that genome.

## Spinach (Spinacia oleracea) isolates

Three spinach-associated genomes were evaluated for redundancy due to
similar strain naming and comparable assembly quality (N50 4.07-4.38 Mb):

- GCA_013347355.2 (Fus187): Skagit County, Washington
- GCA_013347345.2 (Fus254): Washington County, Oregon
- GCA_013347535.2 (Fus167): Arkansas

Each isolate has a distinct NCBI BioSample accession and was collected
from a different U.S. state, indicating three independently sampled
field isolates rather than a single population. All three were retained
to support within-forma-specialis comparison across geographically
distinct populations.

## Strawberry (Fragaria x ananassa) isolates

Five candidate strawberry-associated genomes were identified, each with a
distinct BioSample accession:

- GCA_016164145.2 (MAFF727510): Nara prefecture, Japan
- GCA_016166325.2 (BRIP62122a): Palmview, Queensland, Australia
- GCA_016166095.2 (GL1315): Santa Maria, California, USA
- GCA_016170085.2 (GL1080): Oxnard, California, USA
- GCA_016170095.2 (GL1381): San Luis Obispo, California, USA

Three of the five isolates (GL1315, GL1080, GL1381) originate from
neighboring localities within California's Central Coast strawberry-
growing region, a single connected production area. To avoid over-
representing one regional population, only the highest-contiguity isolate
from this region (GL1381, N50 4.36 Mb) was retained. The Japan and
Australia isolates were retained as independent, geographically distinct
representatives. The final dataset therefore includes three strawberry
genomes spanning three continents.

## Cotton, lettuce, and tomato isolates

The following groups were reviewed for redundancy and retained in full,
as each contains isolates with distinct BioSample accessions from
geographically or lineage-distinct sources:

- Cotton (Gossypium hirsutum): GCA_049306905.1 (VCG 01112, New South
  Wales, Australia) and GCA_049307005.1 (VCG 01111, Queensland,
  Australia) represent two distinct vegetative compatibility groups
  within f. sp. vasinfectum.
- Lettuce (Lactuca sativa): GCA_055944325.1 (California, USA) and
  GCA_055944345.1 (Spain) represent geographically independent isolates.
- Tomato (Solanum lycopersicum): GCA_003315725.1 (Fol4287; race 2, VCG
  0030; Murcia, Spain) and GCA_003977725.1 (isolate D11; race 3; Yolo
  County, California, USA) represent two distinct races within f. sp.
  lycopersici, supporting within-forma-specialis race comparison.

## Forma specialis verification against NCBI taxonomy records

Forma specialis assignments were cross-checked against each assembly's
own submitted NCBI organism name, rather than assumed from host plant
alone. One notable discrepancy was identified: GCA_054643655.1, isolated
from cowpea (Vigna unguiculata), is submitter-designated Fusarium
oxysporum f. sp. phaseoli (the bean-associated lineage) rather than f.
sp. tracheiphilum (the lineage typically associated with cowpea wilt).
This designation is retained as submitted, as it may reflect either a
genuine case of host range extending beyond the lineage's namesake host,
or isolate misclassification at the time of submission; this ambiguity
is noted here rather than resolved by assumption.

For four other genomes (Brassica/conglutinans, banana isolate II5/
cubense, tobacco/nicotianae, potato/tuberosi), the forma specialis shown
in the metadata table reflects published host-association literature
rather than a designation present in the assembly's own NCBI taxonomy
record, and is marked as such in the underlying data.
