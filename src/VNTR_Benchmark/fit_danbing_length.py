"""danbing-tk length predictions for the test individuals.

danbing-tk reports k-mer counts per RPGG locus; its published length estimates
come from a per-locus regression of length on read k-mer dosage, learned from
genomes with both assemblies and short reads.  We follow that design:

  x  = sum of a sample's k-mer counts at the locus / the sample's genome-wide
       MAPQ>=20 depth (from our extraction's per-sample statistics)
  y  = vamos combined length in bp
  per locus, OLS y = a + b x on the training individuals (>= 10 with truth
  and b > 0); otherwise a ratio estimate y = x * mean(y) / mean(x).
  Predictions are made for every sample, and scored only on test individuals.

danbing-tk may widen locus boundaries and merge loci within 700 bp; an RPGG
locus is mapped to our locus only when it overlaps exactly one of ours and that
locus overlaps only it.  Unmapped loci are counted, not scored.

Writes work/danbing/length_predictions.parquet and locus_map.tsv (beside --out).

    python src/VNTR_Benchmark/fit_danbing_length.py
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))
from pipeline.model.assemble import _clean_vamos_col  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402

PIPE = ROOT / "src" / "VNTR_AbsoluteGenotyping" / "work"
DB = ROOT / "src" / "VNTR_Benchmark" / "work" / "danbing"
MIN_FIT = 10


def kmers_per_locus() -> np.ndarray:
    n, cur = [], -1
    with open(DB / "rpgg_290" / "pan.tr.kmers") as fh:
        for line in fh:
            if line[0] == ">":
                n.append(0)
                cur += 1
            elif line.strip():
                n[cur] += 1
    return np.asarray(n)


def locus_sums(path: Path, nk: np.ndarray) -> np.ndarray:
    counts = pl.read_csv(path, has_header=False, new_columns=["c"], schema_overrides={"c": pl.Int64})["c"].to_numpy()
    if counts.size != nk.sum():
        raise ValueError(f"{path.name}: {counts.size} counts vs {nk.sum()} graph k-mers")
    starts = np.concatenate([[0], np.cumsum(nk)[:-1]])
    out = np.zeros(nk.size)
    nz = nk > 0
    out[nz] = np.add.reduceat(counts, starts[nz])
    return out


def map_loci() -> pl.DataFrame:
    # columns 1-3: reference coordinates; then (contig, start, end, strand) per haplotype
    bed = pl.read_csv(DB / "rpgg_290" / "pan.tr.mbe.v2.bed", separator="\t", has_header=False,
                      infer_schema_length=0)
    bed = bed.select(chrom=pl.col(bed.columns[0]), s=pl.col(bed.columns[1]).cast(pl.Int64),
                     e=pl.col(bed.columns[2]).cast(pl.Int64)).with_row_index("dlocus")
    ours = pl.read_csv(DB / "tr_index.tsv", separator="\t")
    by_chrom = defaultdict(list)
    for r in ours.iter_rows(named=True):
        by_chrom[r["chrom"]].append((r["array_start0"], r["array_end0"], r["locus_id"]))
    hits_d, hits_o = defaultdict(list), defaultdict(list)
    for r in bed.iter_rows(named=True):
        for s, e, lid in by_chrom.get(r["chrom"], []):
            if s < r["e"] and e > r["s"]:
                hits_d[r["dlocus"]].append(lid)
                hits_o[lid].append(r["dlocus"])
    rows = [(d, ls[0]) for d, ls in hits_d.items() if len(ls) == 1 and len(hits_o[ls[0]]) == 1]
    m = pl.DataFrame(rows, schema=["dlocus", "locus_id"], orient="row")
    print(f"RPGG loci {bed.height:,}; our loci {ours.height:,}; mapped 1:1 {m.height:,}; "
          f"RPGG loci spanning >1 of ours {sum(len(v) > 1 for v in hits_d.values()):,}; "
          f"our loci with no RPGG locus {ours.height - len(hits_o):,}")
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DB / "length_predictions.parquet",
                    help="predictions parquet (compare_tools.py reads the default path)")
    args = ap.parse_args()
    cfg = load_config(None)
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg)}
    split = dict(pl.read_parquet(PIPE / "model" / "sample_split.parquet").select("sample_id", "fold").iter_rows())
    nk = kmers_per_locus()
    samples = sorted(p.name.split(".")[0] for p in (DB / "align_290").glob("*.done"))
    depth = {s: json.loads((PIPE / "features" / f"{s}.nuisance.json").read_text())["depth_mean"] for s in samples}
    X = np.vstack([locus_sums(DB / "align_290" / f"{s}.tr.kmers", nk) / depth[s] for s in samples])
    print(f"{len(samples)} samples, {nk.size:,} RPGG loci")

    m = map_loci()
    m.write_csv(args.out.with_name("locus_map.tsv"), separator="\t")
    gt = (pl.read_parquet(PIPE / "ground_truth" / "gt_long.parquet", columns=["locus_id", "vamos_sample", "combined_bp", "n_haps_called"])
          .filter(pl.col("n_haps_called") == 2)
          .join(pl.DataFrame({"sample_id": samples, "vamos_sample": [smap[s] for s in samples]}), on="vamos_sample")
          .join(m, on="locus_id"))
    Y = np.full_like(X, np.nan)
    si = {s: i for i, s in enumerate(samples)}
    for r in gt.iter_rows(named=True):
        Y[si[r["sample_id"]], r["dlocus"]] = r["combined_bp"]
    train = np.array([split.get(s) == "train" for s in samples])

    pred = np.full_like(X, np.nan)
    kind = np.full(nk.size, "", dtype=object)
    for d in m["dlocus"].to_list():
        x, y = X[train, d], Y[train, d]
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() >= MIN_FIT and np.ptp(x[ok]) > 0:
            b, a = np.polyfit(x[ok], y[ok], 1)
            if b > 0:
                pred[:, d] = np.clip(a + b * X[:, d], 0, None)
                kind[d] = "ols"
                continue
        if ok.sum() >= 3 and x[ok].mean() > 0:
            pred[:, d] = X[:, d] * y[ok].mean() / x[ok].mean()
            kind[d] = "ratio"
    lid = dict(m.iter_rows())
    out = pl.DataFrame([(lid[d], s, float(pred[i, d]), float(X[i, d]), kind[d])
                        for d in m["dlocus"].to_list() for s, i in si.items() if np.isfinite(pred[i, d])],
                       schema=["locus_id", "sample_id", "danbing_bp", "kmer_dosage", "fit"], orient="row")
    out.write_parquet(args.out)
    print(f"fits: ols {int((kind == 'ols').sum()):,}, ratio {int((kind == 'ratio').sum()):,}; "
          f"wrote {out.height:,} predictions")


if __name__ == "__main__":
    main()
