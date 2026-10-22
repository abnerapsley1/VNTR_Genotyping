#!/bin/bash
# Refit transfer + final models without the 24 externally precomputed population features
# (ANOVA tables, population columns of the mutation-rate table), into work/model_noext.
set -uo pipefail
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
export VNTR_GBM=xgboost
export VNTR_EXCLUDE_FEATURES=$(python -c "
import polars as pl
d=pl.read_csv('docs/feature_inventory.csv')
print(','.join(sorted(set(d.filter(pl.col('origin_note').str.contains('may include'))['feature']))))")
echo "$(date) excluding: $VNTR_EXCLUDE_FEATURES"
rm -rf work/model_noext && cp -a work/model work/model_noext
sed "s|work/model$|work/model_noext|" config/config.yaml > work/model_noext/config.yaml
for step in transfer global evaluate; do
  echo "$(date) $step"
  python -u -m pipeline.cli $step --config work/model_noext/config.yaml || { echo "$(date) FAILED at $step"; exit 1; }
done
echo "$(date) ablation done"
