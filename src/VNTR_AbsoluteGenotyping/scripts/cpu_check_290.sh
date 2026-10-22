#!/bin/bash
# One-time check of the GPU refit: rerun the GBM stages with sklearn on the CPU
# into work/model_cpu_check (starting from the GPU run's upstream outputs).
set -uo pipefail
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
rm -rf work/model_cpu_check && cp -a work/model work/model_cpu_check
sed "s|work/model$|work/model_cpu_check|" config/config.yaml > work/model_cpu_check/config.yaml
for step in transfer global per-allele evaluate; do
  echo "$(date) CPU check: $step"
  VNTR_GBM=sklearn nice -n 10 python -u -m pipeline.cli $step --config work/model_cpu_check/config.yaml \
    || { echo "$(date) CPU check FAILED at $step"; exit 1; }
done
echo "$(date) CPU check done"
