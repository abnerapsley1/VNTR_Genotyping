#!/bin/bash
# Build the danbing-tk RPGG from the 209 training individuals' assemblies with
# bounded disk: genomes are aligned + lifted in batches as their assemblies
# arrive, and each batch's reference-to-assembly PAFs (~3 GB per haplotype) are
# deleted once its liftover is done.  When all genomes are lifted, the joint
# annotation and graph steps run.
set -uo pipefail
cd "$(dirname "$0")/../.."
ROOT=$PWD
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_benchmark
export PATH=$PATH:$ROOT/tools/danbing-tk/bin
ASM=$ROOT/data/Benchmark/assemblies
OUT=$ROOT/src/VNTR_Benchmark/work/danbing/rpgg_290
GEN=$ROOT/src/VNTR_Benchmark/work/danbing/genomes_290.txt
SMK=$ROOT/tools/danbing-tk/pipeline/GoodPanGenomeGraph.snakefile
cp $ROOT/src/VNTR_Benchmark/danbing/goodPanGenomeGraph_290.json $OUT/goodPanGenomeGraph.json
cd $OUT
snakemake -s $SMK --unlock >/dev/null 2>&1 || true   # clear a lock left by a killed run
run_smk() { snakemake -p -s $SMK --cores 32 --resources cores=32 mem=450 --rerun-incomplete --keep-going "$@"; }
while true; do
  todo=(); ready=()
  for g in $(cat $GEN); do
    [ -f $g/lift.foo ] && continue
    todo+=($g)
    [ -f $ASM/$g.0.fa.fai ] && [ -f $ASM/$g.1.fa.fai ] && ready+=($g)
  done
  echo "$(date) lifted: $(( $(wc -l < $GEN) - ${#todo[@]} )) / $(wc -l < $GEN); ready now: ${#ready[@]}"
  if [ ${#todo[@]} -eq 0 ]; then break; fi
  if [ ${#ready[@]} -eq 0 ]; then sleep 600; continue; fi
  batch=("${ready[@]:0:40}")
  targets=$(for g in "${batch[@]}"; do printf '%s/%s/lift.foo ' $OUT $g; done)
  run_smk $targets || echo "$(date) batch had failures; will retry"
  for g in "${batch[@]}"; do
    [ -f $g/lift.foo ] && rm -f r2a/$g.0.r2a.paf r2a/$g.1.r2a.paf
  done
done
echo "$(date) all genomes lifted; running joint annotation + graph"
run_smk && echo "$(date) RPGG_290 built" || echo "$(date) RPGG_290 build FAILED"
