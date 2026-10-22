# Per-haplotype VNTR lengths from short reads — progress report

*Written 2026-10-14/15 while working unattended. All numbers are on held-out people
(validation people for tuning; the 52 test people for the scores below) at repeats
the model was trained on ("known repeats"), unless stated otherwise.*

## The question

Read depth measures the **sum** of a person's two alleles. A total of 40 repeats
could be 20 + 20 or 10 + 30. To split it, extra information is needed: either
knowledge of which allele combinations actually occur, or reads that see one
allele at a time. We tried four sources, alone and together.

"Both alleles exactly right" below means both of the person's allele lengths
(in whole repeat units) match the long-read (vamos) truth. About 20% of
person-sites have two different alleles; those are the hard cases.

## Results on the 52 test people (known repeats, 2.8 million person-sites)

| Method | Both alleles exactly right | Both right where the alleles differ | Both within 1 repeat | "Different alleles" calls: found / correct |
|---|---|---|---|---|
| Equal halves (old per-allele option) | 75.6% | 27.9% | 90.6% | – |
| 1. Allele catalog + allele-gap evidence | 81.7% | 36.8% | 91.7% | 60% / 73% |
| **2. Flanking SNP haplotypes** | **92.6%** | **71.1%** | **96.0%** | **88% / 93%** |
| 3a. Allele catalog + spanning reads | 83.5% | 41.6% | 91.9% | 64% / 77% |
| 3b. SNP haplotypes + spanning reads | 92.4% | 71.0% | 96.1% | 88% / 92% |

Validation people show the same picture (equal halves 74.7%, catalog 80.2%,
SNP haplotypes 91.9%, catalog + reads 82.2%, SNP + reads 91.8%).

**Bottom line:** flanking SNP haplotypes are the method to use. They give both
allele lengths exactly right at 92.6% of person-sites (97.9% at short repeats,
87.4% at longer ones) and correctly phase the alleles. Spanning reads are a
useful backup when phased SNPs are not available (catalog 81.7% → 83.5%; at
short repeats 91.5% → 95.0%), but add nothing on top of SNP haplotypes.

### 1. Allele catalog (`pipeline/model/haplotypes.py`, CLI `haplotypes`)

For each repeat, the panel (long-read people, excluding validation/test people)
gives the alleles that actually occur and their frequencies. Given our predicted
total, the method picks the most likely pair of real alleles that adds up to it.
It also uses a small classifier that predicts how far apart the two alleles are
from fragment-length patterns (spanning read-pair insert sizes, already in our
extracted features). This "allele-gap" evidence added about 1 point.

Outputs per person and repeat: the two allele lengths, the probability that the
two alleles differ, and a 90% range for each allele. Both true alleles fall in
their ranges 93% of the time.

**Repeats the model has never seen** (no catalog allowed): no gain over equal
halves (49.9% both exact either way). The fragment-length signal alone is too
weak to call differing alleles reliably; only 29% of such calls were right.

### 2. Flanking SNP haplotypes (`pipeline/model/haplotypes_snp.py`, CLI `haplotypes-snp --no-reads`)

A repeat allele is usually inherited together with the SNPs around it. All 290
people have phased SNPs from 1000 Genomes (2022 phased panel; chr6 and chr9 from
the 2020 release because EBI's 2022 files fail their own checksums).

- For each repeat, 256 SNPs on each side (within 100 kb) are used.
- Each of a new person's two SNP haplotypes is compared with the 424 training
  haplotypes by how far the identical stretch extends outward from the repeat.
- The 2 best-matching training haplotypes vote on which repeat allele that
  chromosome copy carries.
- Training people's repeat alleles are known, but not which SNP haplotype carries
  which. An EM step works that out first, using only training people.

This replaces the allele-frequency prior in method 1 (all other evidence is the
same), and it also says **which chromosome copy carries which allele** (phase).

### 3. Spanning reads (`pipeline/features/spanning_alleles.py`)

When an allele plus some flanking sequence fits inside one 150 bp read, a read
crossing the whole repeat measures that allele's length directly. That covers
27,920 repeats, about half of all person-sites (allele lengths up to 120 bp).

- Reads are internally consistent: where a person has two equal alleles, about
  100% of reads agree on one length.
- Our repeat boundaries differ from vamos's by a small fixed amount at each repeat
  (0 bp at most; sometimes ±1–4 bp). This offset is learned per repeat from 15
  training people whose CRAMs are on disk.
- Extracted for the 52 test people (kept CRAMs) and the 26 validation people
  (CRAMs downloaded from AWS, MD5-checked, then deleted). Reading the AWS copies
  in place silently lost regions when connections dropped, so a person is now
  rejected unless at least 95% of repeats have reads (complete files reach about
  98%).
- Read evidence is weighted against the other evidence. The weight was tuned on
  validation people.

Results at the short repeats where spanning reads exist (1.39 million test
person-sites; 12% have two different alleles):

| Method | Both exactly right | Where the alleles differ | "Different" calls: found / correct |
|---|---|---|---|
| Catalog | 91.5% | 56.5% | 60% / 71% |
| Catalog + spanning reads | 95.0% | 72.8% | 75% / 87% |
| SNP haplotypes | **97.9%** | **87.7%** | **90% / 95%** |
| SNP haplotypes + spanning reads | 97.6% | 87.2% | 89% / 93% |

At the remaining (longer) repeats: SNP haplotypes 87.4% (64.4% where the
alleles differ) vs catalog 72.1% (28.9%).

### 4. Family checks (`scripts/family_check.py`)

The test set has one complete trio (SH032: HG00512 father, HG00513 mother,
HG00514 daughter) and a father–son duo (GB40). None of them informed any model.

| Method | Trio: child's alleles inheritable from parents | …where the child's alleles differ | Duo: share an allele |
|---|---|---|---|
| Long-read truth (vamos) | 98.7% | 95.3% | 99.9% |
| SNP haplotypes | 97.6% | 90.1% | 98.4% |
| SNP haplotypes + spanning reads | 97.5% | 89.7% | 98.4% |
| Allele catalog + spanning reads | 96.4% | 86.9% | 98.6% |
| Allele catalog | 96.1% | 85.9% | 98.5% |
| Equal halves | 86.3% | 57.4% | 90.8% |

**Phase check.** The child's SNP haplotypes are pedigree-phased, so each one comes
from one parent along a whole chromosome. Where only one parent carries the
child's longer allele, the method's phase call agrees with that parent 99.0% of
the time (3,458–3,484 informative repeats depending on the variant; confident
calls only). On all 22 chromosomes
the child's first haplotype is the father's, which is the 1000 Genomes convention.

### Caveats

- The SNP method needs phased SNP haplotypes for the person. 1000 Genomes people
  have them. For a new sample, they would come from statistically phasing its own
  SNP calls against a reference panel. That is standard, but a little less
  accurate than the pedigree-aided phasing used here.
- The catalog and SNP methods only know alleles seen in the long-read panel.
  Alleles outside the panel fall back to equal halves.
- New repeats (outside the panel) get no per-haplotype improvement yet.

---

# Our tool vs danbing-tk: statistical tests

*Script: `src/VNTR_Benchmark/stats_ours_vs_danbing.py` and `stats_supplement.py`;
tables in `src/VNTR_Benchmark/work/compare/stats_ours_vs_danbing.csv` and
`stats_supplement.csv`.*

Both tools are compared on exactly the same calls: 2.63 million person-sites (52
test people × 51,695 repeats) where both made a call. Calls from the same person
or the same repeat are not independent, so the tests use **people** or
**repeats** as units, never individual calls:

- **Per person:** each of the 52 people gets one score per tool; Wilcoxon
  signed-rank test on the 52 paired differences.
- **Bootstrap:** 95% confidence intervals for the difference, resampling whole
  people, and separately whole repeats.
- **Per repeat:** for tracking person-to-person differences, a correlation per
  repeat (42,639 repeats where people differ); Wilcoxon on the paired values.
- **Rare alleles:** McNemar exact test on the paired right/wrong calls.
- **Multiple testing:** p-values Holm-corrected across all tests.

| Measure | Our tool | danbing-tk | Difference (95% CI, resampling people) | People where ours is better | p (Holm) |
|---|---|---|---|---|---|
| Median error (% of length) | 0.80% | 1.77% | −0.97 points (−0.99 to −0.95) | 52 / 52 | 2×10⁻⁸ |
| Calls within 5% | 83.1% | 68.9% | +14.2 points (14.1 to 14.5) | 52 / 52 | 2×10⁻⁸ |
| Mean absolute error | 22.4 bp | 38.9 bp | −16.5 bp (−17.2 to −15.8) | 52 / 52 | 2×10⁻⁸ |
| Tracks person-to-person differences (r) | 0.90 | 0.79 | +0.11 (0.10 to 0.12) | 33,742 vs 8,897 repeats | < 10⁻³⁰⁰ |
| Rare alleles, right direction | 72.9% | 63.7% | +9.2 points (7.1 to 10.8) | 370 vs 156 calls | 8×10⁻²⁰ |

(p = 2×10⁻⁸ is about the smallest possible here: all 52 of 52 people favour our tool.)

- **Every group:** ours is significantly better on accuracy and tracking (p < 0.001
  after correction) in every core-spread group, site type, array-length group,
  and both inside and outside segmental duplications.
- **Rare alleles at the most variable repeats (core spread ≥1,000 bp):** danbing-tk
  had looked better here (83% vs 53%), but that compared different site sets. On
  calls both tools made there are only 12 rare-allele cases, and the two tools got
  exactly the same ones right, so there is no evidence of a difference.
- **No significant difference** (after correction) for rare alleles at stable
  sites, at sites with core spread 300–999 bp, and in segmental duplications.
  These groups have few rare-allele cases.
