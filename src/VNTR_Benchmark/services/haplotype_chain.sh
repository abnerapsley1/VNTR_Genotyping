#!/bin/bash
# 1) SNP files complete -> SNP-haplotype model (no reads) + family check
# 2) spanning-read extraction complete -> full model (SNP + reads, with ablations) + family check
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
L=work/logs/haplotype_chain.log
until grep -q "all chromosomes done" /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_Benchmark/work/logs/fetch_snps.log; do sleep 120; done
echo "$(date) SNPs ready: $(ls /data/GitHub_Repos/VNTR_Genotyping/data/1000Genomes/phased_snps/*.csi | wc -l) chromosome files" >> $L
python -u -m pipeline.cli haplotypes-snp --no-reads >> $L 2>&1 && cp work/model/haplotype_snp_eval.json work/model/haplotype_snp_eval.noreads.json \
  && cp work/model/haplotype_snp_predictions.parquet work/model/haplotype_snp_predictions.noreads.parquet
python -u scripts/family_check.py >> $L 2>&1
until grep -q "spanning extraction done" work/logs/spanning.log; do sleep 120; done
echo "$(date) spanning reads ready: $(ls work/features/spanning/*.parquet | wc -l) people" >> $L
python -u -m pipeline.cli haplotypes-snp >> $L 2>&1
python -u scripts/family_check.py >> $L 2>&1
echo "$(date) haplotype chain done" >> $L
