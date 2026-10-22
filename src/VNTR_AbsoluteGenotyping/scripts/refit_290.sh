#!/bin/bash
# After feature extraction finishes for all 290 samples:
#   1. keep the 70-sample model outputs in work/model_70
#   2. refit the full model chain on the GPU (VNTR_GBM=xgboost) into work/model
#   3. one-time check: rerun the GBM stages on the CPU (sklearn) into
#      work/model_cpu_check, at low priority, so the two evaluations can be compared
set -uo pipefail
cd /data/GitHub_Repos/VNTR_Genotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
cd src/VNTR_AbsoluteGenotyping
while systemctl --user is-active --quiet vntr-extract.service; do sleep 600; done
echo "$(date) extraction finished; retrying any samples that failed"
python -u -m pipeline.cli extract     # resumes: completed samples are skipped
n=$(ls work/features/*.parquet | wc -l)
echo "$(date) $n feature tables"
[ "$n" -ge 290 ] || { echo "$(date) only $n of 290 samples extracted; not refitting"; exit 1; }
[ -d work/model_70 ] || cp -a work/model work/model_70
echo "$(date) GPU refit"
VNTR_GBM=xgboost python -u -m pipeline.cli all-model || { echo "$(date) GPU refit FAILED"; exit 1; }
echo "$(date) GPU refit done"
rm -rf work/model_cpu_check && cp -a work/model work/model_cpu_check
sed "s|work/model$|work/model_cpu_check|" config/config.yaml > work/model_cpu_check/config.yaml
for step in transfer global per-allele evaluate; do
  VNTR_GBM=sklearn nice -n 10 python -u -m pipeline.cli $step --config work/model_cpu_check/config.yaml \
    || { echo "$(date) CPU check FAILED at $step"; exit 1; }
done
echo "$(date) CPU check done"
