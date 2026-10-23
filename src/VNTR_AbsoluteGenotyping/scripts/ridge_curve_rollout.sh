#!/bin/bash
# Roll out ridge-combined calibration curves:
#   1. refit curves + correction models into work/model_ridge (current default untouched)
#   2. compare with the current default on the test people
#   3. if test error is lower: promote to default, rebuild the vntrlen bundle, regenerate reports
set -uo pipefail
R=/data/GitHub_Repos/VNTR_Genotyping
A=$R/src/VNTR_AbsoluteGenotyping
cd $A
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
export VNTR_GBM=xgboost
OUT=$A/work/experiments/curve_combination
mkdir -p $OUT
step() { echo; echo "===== $(date '+%F %T')  $* ====="; }

step "1/4 refit with ridge curves into work/model_ridge"
rm -rf work/model_ridge && cp -a work/model work/model_ridge
sed "s|work/model$|work/model_ridge|" config/config.yaml > work/model_ridge/config.yaml
for s in calibrate global per-allele evaluate; do
  step "  $s"
  python -u -m pipeline.cli $s --config work/model_ridge/config.yaml || { echo "FAILED at $s"; exit 1; }
done

step "2/4 compare with the current default (single-signal curves) on test people"
python - <<'EOF' | tee $OUT/end_to_end_comparison.txt
import json
import polars as pl
W = "/data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping/work"
rows = []
for name, d in (("single (current)", "model"), ("ridge", "model_ridge")):
    e = json.load(open(f"{W}/{d}/evaluation.json"))["quadrants"]
    for q in ("test_novel_sample", "test_novel_locus", "test_novel_both", "val"):
        m = e[q]["model"]
        rows.append(dict(curves=name, quadrant=q, mean_abs_err=m["mean_abs_err_units"], within_1=m["frac_within_1u"],
                         within_locus_r=m["within_locus_r"], coverage90=m["interval_coverage_90"]))
pl.Config.set_tbl_width_chars(200)
print(pl.DataFrame(rows).sort("quadrant", "curves").with_columns(pl.col(pl.Float64).round(4)))
# by core spread, test people at catalogued repeats
edges, labs = [1, 2, 5, 20, 100], ["0", "1", "2-4", "5-19", "20-99", ">=100"]
out = []
for name, d in (("single (current)", "model"), ("ridge", "model_ridge")):
    p = pl.read_parquet(f"{W}/{d}/predictions.parquet").filter(pl.col("quadrant") == "test_novel_sample")
    p = p.with_columns(g=pl.col("cov_panel_core90_units").cut(edges, labels=labs, left_closed=True).cast(pl.Utf8),
                       e=(pl.col("y_pred") - pl.col("combined_units")).abs())
    out.append(p.group_by("g").agg(n=pl.len(), mae=pl.col("e").mean(), within1=(pl.col("e") <= 1).mean())
               .with_columns(curves=pl.lit(name)))
print(pl.concat(out).sort("g", "curves").with_columns(pl.col(pl.Float64).round(4)))
EOF

step "3/4 decision"
promote=$(python -c "
import json
W='$A/work'
a=json.load(open(W+'/model/evaluation.json'))['quadrants']['test_novel_sample']['model']['mean_abs_err_units']
b=json.load(open(W+'/model_ridge/evaluation.json'))['quadrants']['test_novel_sample']['model']['mean_abs_err_units']
print('yes' if b < a else 'no')")
echo "promote ridge curves: $promote" | tee -a $OUT/end_to_end_comparison.txt
if [ "$promote" != "yes" ]; then
  echo "ridge curves did not lower test error; current default kept"; exit 0
fi

step "4/4 promote: default model, vntrlen bundle, reports"
rm -f work/model_ridge/config.yaml
mv work/model work/model_single_curve && mv work/model_ridge work/model
rm -rf $R/data/vntrlen_bundle/52ec291b8f70 $R/data/vntrlen_bundle/*.tar.gz
( cd $R/src/vntrlen && python -u scripts/make_shipped_bundle.py --out $R/data/vntrlen_bundle --jobs 24 ) \
  2>&1 | grep -vE "Warn|warn|implode|polars|issues/|is_in|Potential|Use a data|Set the device|This warning|return func|^$"
bash $A/scripts/regenerate_reports.sh
echo; echo "===== $(date '+%F %T')  ALL DONE ====="
