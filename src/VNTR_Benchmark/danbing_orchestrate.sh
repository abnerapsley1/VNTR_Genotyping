#!/bin/bash
# danbing-tk end to end, unattended:
#   1. wait for all 96 training-assembly haplotypes (.fa + .fai)
#   2. build the RPGG with danbing-tk's GoodPanGenomeGraph pipeline (local, 32 cores)
#   3. genotype every benchmark CRAM (15 test first, then 48 training) as it arrives
set -euo pipefail
cd "$(dirname "$0")/../.."
ROOT=$PWD
source /data/miniconda3/etc/profile.d/conda.sh
conda activate vntr_benchmark
export PATH=$PATH:$ROOT/tools/danbing-tk/bin
ASM=data/Benchmark/assemblies
RPGG=src/VNTR_Benchmark/work/danbing/rpgg

echo "$(date) waiting for 96 assembly haplotypes"
until [ "$(ls $ASM/*.fa.fai 2>/dev/null | wc -l)" -ge 96 ] && ! pgrep -f "^python -u src/VNTR_Benchmark/fetch_assemblies" >/dev/null; do sleep 300; done
echo "$(date) assemblies ready; building RPGG"

if [ ! -s $RPGG/pan.kmerDBi.vv ]; then
  cd $RPGG
  cp $ROOT/src/VNTR_Benchmark/danbing/goodPanGenomeGraph.json .
  snakemake -p -s $ROOT/tools/danbing-tk/pipeline/GoodPanGenomeGraph.snakefile \
      --cores 32 --resources cores=32 mem=450 --rerun-incomplete --keep-going
  cd $ROOT
fi
echo "$(date) RPGG built"

conda activate vntr_analysis
TEST="HG00513 HG00621 HG00733 HG01071 HG01106 HG01109 HG01175 HG01505 HG02011 HG02148 HG02717 HG02818 HG03732 NA18906 NA19240"
TRAIN=$(cat src/VNTR_Benchmark/work/danbing/genomes.txt)
# two samples at a time: one's name-collation (~2 cores) overlaps the other's alignment
while true; do
  ready=(); left=0
  for s in $TEST $TRAIN; do
    [ -f src/VNTR_Benchmark/work/danbing/align/$s.done ] && continue
    left=$((left+1))
    [ -f data/Benchmark/cram/$s.cram ] && ready+=("$s")
  done
  [ $left -eq 0 ] && { echo "$(date) all danbing-tk alignments done"; exit 0; }
  if [ ${#ready[@]} -gt 0 ]; then
    printf '%s\n' "${ready[@]}" | xargs -P 2 -I{} python -u src/VNTR_Benchmark/run_danbing_align.py {} --threads 16 || true
  else
    sleep 300
  fi
done
