#!/bin/bash
# redo validation people whose spanning-read extraction was incomplete, then the full haplotype model
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
L=work/logs/haplotype_final.log
python -u -m pipeline.features.spanning_alleles --fold val --jobs 28 >> $L 2>&1
n=$(ls work/features/spanning/*.parquet | wc -l)
echo "$(date) spanning files: $n" >> $L
[ "$n" -ge 93 ] || { echo "$(date) incomplete spanning extraction; stopping" >> $L; exit 1; }
python -u -m pipeline.cli haplotypes-snp >> $L 2>&1 || { echo "$(date) haplotypes-snp FAILED" >> $L; exit 1; }
python -u scripts/family_check.py >> $L 2>&1
echo "$(date) haplotype final done" >> $L
