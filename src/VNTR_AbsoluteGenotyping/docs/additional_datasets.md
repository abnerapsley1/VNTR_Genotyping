# Additional public data with both phased long-read assemblies and short reads

Search date: 2026-10-09. Current panel: 70 individuals (36 HPRC + 34 HGSVC2 short-read lists).

**Bottom line:** 220 more individuals are usable right away, with no new truth data needed. They are already in our vamos call set (`vamosExpanded_v3.0_GRCh38_oriMotifs`) and have public 1000 Genomes 30× Illumina CRAMs. That raises the panel from 70 to 290. The other cohorts below either lack public short reads or would need vamos run on their assemblies first.

The 220 individuals, with population, sex, ENA run, CRAM URL, MD5 and relatives already in the panel, are listed in [`additional_samples_vamos_1kgp.tsv`](additional_samples_vamos_1kgp.tsv). The counts in section A were computed directly from the vamos VCF header and the two 1000 Genomes high-coverage sequence indexes (3,202 samples).

## Summary

| Dataset | Individuals with phased assemblies | Usable (assembly + public short reads) | Not yet in our panel | Long-read data | Short-read data | In our vamos truth? | Access |
|---|---|---|---|---|---|---|---|
| **HPRC Release 2** | 232 | 223 (1kGP 30×) + up to 9 more via HPRC/GIAB Illumina | **184** | PacBio HiFi + ONT ultra-long; hifiasm, phased by trio Illumina or Hi-C | NYGC 1kGP 30× Illumina NovaSeq 6000, 150 bp PE | Yes (232 columns) | Open |
| **HGSVC phase 3** | 65 | 59 in vamos with 1kGP 30× | **31** | PacBio HiFi + ONT ultra-long; Verkko, phased by Strand-seq | NYGC 1kGP 30× | Yes (60 columns) | Open |
| **JaSaPaGe (Japanese + Saudi; "JSA" in vamos)** | 14 in vamos (10 Japanese in paper + Saudi) | 5 (Japanese, in 1kGP JPT) | **5** | PacBio HiFi + ONT + Hi-C | NYGC 1kGP 30× for the 5 JPT; Saudi short reads not found publicly | Yes (14 columns) | 1kGP open; Saudi not public |
| **Chinese Pangenome Consortium phase 1** | 58 | 0 public | 0 | PacBio HiFi (30.7× avg) or ONT, linked reads, Hi-C | Illumina generated, but controlled access | Yes (57 columns) | GSA-Human HRA002344, controlled |
| **Arab Pangenome Reference** | 53 | 0 confirmed public | 0 | PacBio HiFi 35×, ONT ultra-long 54×, Hi-C | High-coverage Illumina generated; public access not confirmed | Yes (50 columns, "APR") | Not verified; likely controlled |
| **Platinum Pedigree (CEPH 1463)** | 28-member, 4-generation family; assemblies for a subset (count not stated in repo) | 28 have Illumina; NA12878 (used), NA12891, NA12892 also in 1kGP | Up to ~27, all relatives | PacBio HiFi, ONT ultra-long, Strand-seq, Element | Illumina + Element, open on AWS | No (would need vamos run on assemblies); has its own tandem-repeat truth set | Open (AWS) |
| **1000 Chinese Pangenome** | 1,116 (2,232 haplotypes) | Short and long reads generated; public access not verified | Unknown | Long reads + short reads + Hi-C (platform not stated in news coverage) | Generated | No | Not verified; Chinese genomic data is usually controlled |
| **Korean Pangenome Reference (K-PanRef)** | 13 diploid + KOREF1 | Unknown | Unknown | PacBio HiFi (national project) | Unknown | No | Not verified |
| **1kGP ONT (Schloissnig et al. 2025)** | **0**: no phased assemblies | 1,019 with 1kGP 30× | Silver-standard only | ONT, 16.9× median, read N50 20.3 kb | NYGC 1kGP 30× | No, but the paper ran vamos on reads (370,424 VNTRs) | Open (1KG_ONT_VIENNA) |

## A. Ready to add now (vamos truth already exists, open 30× CRAMs)

| vamos cohort | vamos columns | With 1kGP 30× | Already in our 70 | **Additional** |
|---|---|---|---|---|
| HPRC phase 2 | 232 | 223 | 39 | **184** |
| HGSVC phase 3 | 60 | 59 | 28 | **31** |
| JSA (Japanese part) | 14 | 5 | 0 | **5** |
| HPRC phase 1 / HGSVC phase 2 | 3 | 3 | 3 | 0 |
| APR, CPC | 107 | 0 | 0 | 0 |
| **Total** | **416** | **290** | **70** | **220** |

**Ancestry of the 220** (1kGP superpopulation): AFR 63, EAS 54, SAS 39, EUR 32, AMR 32. The combined 290 would be AFR 91, EAS 65, AMR 51, SAS 45, EUR 38. That is much more balanced than most references, and Africa-rich, which suits VNTRs: African genomes carry the most length diversity.

**Short-read uniformity.** All 290 come from the same NYGC 30× NovaSeq 6000 run (BWA-MEM, GRCh38 analysis set with decoys and HLA). Sequencing differences won't confound the panel, but the calibration also won't learn other protocols.

**Things to know before adding them:**
- **Relatives.** The 1kGP set includes trios. 11 parent–child pairs fall within the 290 (e.g. HG00731/HG00732/HG00733, NA19238/NA19239/NA19240, HG01890/HG01891), and 7 of the 220 new individuals have a parent or child in the panel (see the `relatives_in_panel` column). Sample folds should be assigned **by family**, so a child and parent never split between training and test. Relatives share alleles and would inflate test scores.
- **Truth comes from different assembly releases.** HPRC phase 2 (hifiasm trio/Hi-C) and HGSVC phase 3 (Verkko) differ in assembler and phasing. The current panel already mixes both, so this doesn't change the design.
- **Nine more individuals** in vamos (HPRC/HGSVC: HG002, HG005, HG01123, HG02109, HG02486, HG02559, HG03471, HG06807, NA21309/NA21487) are outside the 3,202-sample 1kGP set. HPRC Release 2 lists Illumina WGS for over 200 of its samples, and HG002/HG005 have GIAB Illumina. They are probably usable after checking those sources, but that wasn't verified here. Their reads would come from a different run than the NYGC data.
- **Cost.** About 220 × 15–30 GB ≈ 4–5 TB of CRAM downloads. Our download mode deletes each CRAM after extraction, so peak disk use stays low, but downloads take about 2–3 days at current speeds. Extraction plus model refit adds about 1–2 days of compute. Every CRAM should be MD5-verified against these indexes (the benchmark found 8 of 51 downloads corrupted).

## B. Would need extra work or access

- **Chinese Pangenome Consortium phase 1.** 58 individuals from 36 minority ethnic groups. Raw reads and assemblies are under controlled access (GSA-Human HRA002344). Chinese human-genetic-resource rules make export approval slow and uncertain.
- **Arab Pangenome Reference.** 53 individuals from 8 Arab countries; HiFi, ONT ultra-long, Hi-C and high-coverage Illumina were all generated. I could not confirm public access to the Illumina data; the paper's data-availability statement needs checking. It's the only Middle-Eastern/North-African source, so worth a request if available.
- **JaSaPaGe Saudi samples.** 9 Saudi assemblies are in vamos ("ksa" columns); no public short reads were found.
- **Platinum Pedigree (CEPH 1463).** Fully open on AWS, with Illumina, Element, HiFi and ONT data and a pedigree-validated tandem-repeat truth set. All members are relatives of NA12878, so it adds little population diversity. It would, however, be a strong independent check: Mendelian consistency, a second sequencing chemistry (Element), and truth from a different method. vamos would need to be run on its assemblies.
- **1000 Chinese Pangenome (2026) and Korean pangenome.** Large and recent, but public raw-data access isn't confirmed and they aren't in vamos. Worth monitoring.

## C. A different kind of resource: 1kGP ONT read-based VNTR calls

Schloissnig et al. (2025) sequenced 1,019 1kGP individuals with ONT (16.9× median) and genotyped 370,424 VNTRs with vamos **from reads, not assemblies**. All 1,019 have NYGC 30× short reads. These calls are noisier than assembly-based truth, but could serve as:
- a **large secondary test set**, especially for population-level checks
- a **source of rare-allele examples** for weak supervision

It shouldn't replace assembly-based truth for headline accuracy numbers.

## Recommendation

1. **Add the 220 now.** It quadruples the panel with identical data types and no new truth work. Assign folds by family, and MD5-check every CRAM.
2. **Check the 9 non-1kGP vamos individuals** for HPRC or GIAB Illumina data (a small, cheap addition).
3. **Request the Arab Pangenome Illumina data** if it is available on request.
4. **Use Platinum Pedigree as an external validation set** rather than training data.
5. **Consider the 1kGP ONT vamos calls** as a large silver-standard test set.

## Sources

- HPRC Release 2: [humanpangenome.org/hprc-data-release-2](https://humanpangenome.org/hprc-data-release-2/), [IGSR HPRC2 collection](https://www.internationalgenome.org/data-portal/data-collections/hprc2), assembly index [`assemblies_release2_v1.0.index.csv`](https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables)
- HGSVC3: Logsdon et al. 2025, [Complex genetic variation in nearly complete human genomes](https://www.nature.com/articles/s41586-025-09140-6), Nature; [IGSR HGSVC3 collection](https://www.internationalgenome.org/data-portal/data-collections/hgsvc3.html)
- 1kGP high-coverage short reads: [Byrska-Bishop et al., expanded 1000 Genomes cohort including 602 trios](https://www.biorxiv.org/content/10.1101/2021.02.06.430068.full.pdf); sequence indexes `1000G_2504_high_coverage.sequence.index` and `1000G_698_related_high_coverage.sequence.index` (IGSR FTP)
- Chinese Pangenome Consortium: Gao et al. 2023, [A pangenome reference of 36 Chinese populations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10322713/), Nature; data GSA-Human HRA002344
- Arab Pangenome: [Nature Communications 2025, doi:10.1038/s41467-025-61645-w](https://www.genengnews.com/topics/omics/first-arab-human-pangenome-sequenced); [Bio-IT World summary](https://www.bio-itworld.com/news/2025/07/29/arab-pangenome-reference-published)
- JaSaPaGe: [Phased genome assemblies and pangenome graphs of human populations of Japan and Saudi Arabia](https://borg.kaust.edu.sa/node/8471)
- Platinum Pedigree: [Platinum Pedigree Datasets (GitHub)](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets), [AWS Open Data](https://registry.opendata.aws/platinum-pedigree/), [Porubsky et al. bioRxiv](https://www.biorxiv.org/content/10.1101/2024.08.05.606142.full.pdf)
- 1000 Chinese Pangenome: [news summary](https://www.azolifesciences.com/news/20260416/1116-Genomes-Mapped-to-Reveal-Hidden-Human-DNA-Variation.aspx), Nature 2026, doi:10.1038/s41586-026-10315-y
- Korean pangenome: [K-PanRef / KOREF1 preprint](https://www.biorxiv.org/content/10.1101/2025.11.17.688257.full.pdf)
- 1kGP ONT: Schloissnig et al. 2025, [Structural variation in 1,019 diverse humans based on long-read sequencing](https://pmc.ncbi.nlm.nih.gov/articles/PMC12350158/), Nature
- vamos: [Ren et al. 2023, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03010-y)
