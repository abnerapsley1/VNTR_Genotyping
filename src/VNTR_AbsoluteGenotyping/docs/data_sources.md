# Data sources: long-read truth and per-repeat annotation tables

All four files in `data/1000Genomes/` below come from **TRCompDB** (Bida Gu and
Mark Chaisson, Chaisson lab, University of Southern California), the population
tandem-repeat reference built with vamos on 416 diploid long-read assemblies
(832 haplotypes).

| File (ours) | Content | Public release |
|---|---|---|
| `vamosExpanded_v3.0_GRCh38_oriMotifs.vcf` / `.tsv` | vamos genotypes (training truth) and motif catalog | v3.1 files in Zenodo 10.5281/zenodo.21961073 |
| `mutRate_GRCh38_masterFile.tsv` | per-repeat mutation rate + locus features; column definitions in `mutRate_GRCh38_masterFile.columns.txt` | same record (`mutRate_GRCh38_masterFile.tsv.gz`, `.txt`) |
| `anova.annoLen.tsv` | per-repeat ANOVA of allele length across continental groups (EUR, AMR, AFR, SAS, EAS, MDE) | same record |
| `anova.ed2Major.tsv` | per-repeat ANOVA of edit distance to the major allele (allele composition) | same record |

**Version note.** Our copies were received 2026-03-05 and are labelled v3.0; the
public Zenodo record (2026-08-16, CC-BY-4.0) is v3.1, which removed low-motif-purity
loci and has more ANOVA rows (1,228,532 vs our 1,136,502). Results here use v3.0.

**How the tables were made (from the preprint and the column file).**
- Mutation rates: SNP-based local genealogies (ancestral recombination graphs built
  in 2 Mb windows) give pairwise branch lengths between haplotypes; the rate is total
  pairwise edit distance of the repeat alleles divided by total pairwise branch length
  (per generation), optionally normalised by length.  Trio children and poorly covered
  genomes are excluded.  `allGC` uses the 94 HPRC phase 1 haplotypes.
- ANOVAs: per repeat, allele length (or edit distance to the major allele) compared
  across continental groups; MDE = Middle Eastern genomes in the panel.
- Both use the full panel, which includes our validation and test people. For that
  reason the 24 population-summary columns are excluded from the model
  (`EXTERNAL_POPULATION_FEATURES` in `pipeline/model/assemble.py`).

**Citations.**
- Gu B, Peng D, LaFlamme CW, Bennett MF, Bahlo M, Weisburd B, Mefford H, HGSVC, HPRC,
  Chaisson M. *A population-scale map of human tandem repeat composition and mutation
  dynamics from long-read assemblies.* bioRxiv 2024.08.07.607105 (v2, 2026-01-24).
  (v1 title: *TRCompDB: A reference of human tandem repeat sequence and composition
  variation from long-read assemblies*.)
- Gu B, Chaisson M. *TRcompDB v2.1, a global reference of human tandem repeat
  composition and mutation rate from long-read assemblies.* Zenodo (2026).
  doi:10.5281/zenodo.21961073
- Ren J, Gu B, Chaisson MJP. *vamos: variable-number tandem repeats annotation using
  efficient motif sets.* Genome Biology (2023).
