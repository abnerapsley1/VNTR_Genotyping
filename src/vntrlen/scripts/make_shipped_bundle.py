"""Build the shipped GRCh38 model bundle from the research pipeline's outputs.

Training data: the 212 training people (the 26 validation people calibrate the 90%
intervals; the 52 test people are never used).  Every one of the 60,504 repeats gets
calibration curves fitted on the training people, so all repeats go through the same
(known-repeat) pipeline.  Per-repeat performance comes from the 52 test people, scored
with the package's own code (also a check that it reproduces the research predictions).

    python scripts/make_shipped_bundle.py --out /path/to/bundle_dir [--jobs 24]
-> <out>/<fingerprint>/shipped_GRCh38_v1/ and <out>/vntrlen-shipped-GRCh38-v1.tar.gz
"""
import argparse
import json
import sys
import tarfile
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from joblib import Parallel, delayed

PKG = Path(__file__).resolve().parents[1]
RESEARCH = PKG.parent / "VNTR_AbsoluteGenotyping"
sys.path.insert(0, str(PKG))
sys.path.insert(0, str(RESEARCH))            # to unpickle the research models
from vntrlen import checks, curves, gbm, report, stack, store  # noqa: E402
from vntrlen.genotype import add_curve_estimates, tiers  # noqa: E402
from pipeline.model.assemble import _clean_vamos_col  # noqa: E402
from pipeline.model.panel_stats import PANEL_COLUMNS, panel_stats  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402

SET = "shipped_GRCh38_v1"
W = RESEARCH / "work"
MIN_N = 8


def _fit_chunk(items):
    return [(lid, curves.dumps(curves.fit_locus(sig, y, MIN_N))) for lid, sig, y in items]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=24)
    args = ap.parse_args()
    cfg = load_config(None)
    fasta = Path("/data/GitHub_Repos/VNTR_Genotyping/data/Benchmark/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    fp = checks.reference_fingerprint(fasta)
    d = args.out / fp / SET
    (d / "stack").mkdir(parents=True, exist_ok=True)

    # ---- repeats ----
    tg = pl.read_parquet(W / "ground_truth" / "loci_targets.parquet")
    loci = tg.select("locus_id", "chrom", start=pl.col("array_start0"), end=pl.col("array_end0"),
                     period=pl.col("period").cast(pl.Int64), motifs=pl.col("motifs").list.join(","),
                     name=pl.lit(None, dtype=pl.Utf8))
    print(f"{loci.height:,} repeats")

    # ---- annotation + panel statistics (panel without validation/test people, as in training) ----
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg) if s.get("vamos_column")}
    sp = pl.read_parquet(W / "model" / "sample_split.parquet")
    held = {smap[s] for s in sp.filter(pl.col("fold") != "train")["sample_id"] if s in smap}
    cov = pl.read_parquet(W / "covariates" / "locus_covariates.parquet")
    drop = ["chrom", "start", "end", "gene", "locus_key", "match_type"] + PANEL_COLUMNS
    cov = cov.drop([c for c in drop if c in cov.columns])
    cov = cov.join(panel_stats(pl.read_parquet(W / "ground_truth" / "gt_long.parquet"), held), on="locus_id",
                   how="left")
    cov = cov.rename({c: f"cov_{c}" for c in cov.columns if c != "locus_id"})
    cov = cov.select(["locus_id"] + [c for c, t in zip(cov.columns, cov.dtypes) if c != "locus_id" and t.is_numeric()])
    trf = pl.read_parquet(W / "model" / "transfer_predictions.parquet").select("locus_id", "transfer_ood") \
        .unique("locus_id")

    # ---- curves for every repeat, fitted on the training people ----
    sig_cols = curves.SIGNAL_SHORTLIST
    tbl = pl.read_parquet(W / "model" / "model_table.parquet",
                          columns=["locus_id", "sample_id", "sample_fold", "combined_units"] + sig_cols)
    tr = tbl.filter(pl.col("sample_fold") == "train")
    items = [(lid, {c: g[c].cast(pl.Float64).to_numpy() for c in sig_cols}, g["combined_units"].to_numpy())
             for (lid,), g in tr.group_by(["locus_id"])]
    print(f"fitting curves for {len(items):,} repeats on {tr['sample_id'].n_unique()} training people")
    chunks = [items[i::args.jobs * 8] for i in range(args.jobs * 8)]
    fitted = dict(x for part in Parallel(n_jobs=args.jobs)(delayed(_fit_chunk)(c) for c in chunks) for x in part)
    n_tr = dict(tr.group_by("locus_id").len().iter_rows())

    def sig_of(cj):
        c = curves.loads(cj)
        return (c["pe"] or c["se"] or [None])[0]

    bpu = cov.select("locus_id", bp_per_unit=pl.col("cov_mean_combined_bp") / pl.col("cov_mean_combined_units"))
    mlen = loci.select("locus_id", m=pl.col("motifs").str.split(",").list.eval(pl.element().str.len_chars())
                       .list.mean())
    bpu = bpu.join(mlen, on="locus_id").select("locus_id", bp_per_unit=pl.coalesce(
        pl.when(pl.col("bp_per_unit").is_finite()).then(pl.col("bp_per_unit")), pl.col("m")))
    L = (loci.join(cov, on="locus_id", how="left").join(trf, on="locus_id", how="left")
         .join(bpu, on="locus_id", how="left")
         .with_columns(curves=pl.col("locus_id").replace_strict(fitted, default=None),
                       n_train_people=pl.col("locus_id").replace_strict(n_tr, default=0).cast(pl.Float64)))
    L = L.with_columns(curve_signal=pl.col("curves").map_elements(sig_of, return_dtype=pl.Utf8))

    # ---- stack: export the research models to portable boosters ----
    g = joblib.load(W / "model" / "global_model.joblib")["known"]
    t = joblib.load(W / "model" / "transfer_model.joblib")
    for key, fname in (("point", "point.json"), ("q_lo", "q_lo.json"), ("q_hi", "q_hi.json")):
        g[key].model.get_booster().save_model(str(d / "stack" / fname))
    t["model"].model.get_booster().save_model(str(d / "stack" / "transfer.json"))
    ood_cols = [c for c in t["features"] if c.startswith("cov_")]
    train_loci = pl.read_parquet(W / "model" / "locus_split.parquet").filter(pl.col("fold") == "train")["locus_id"]
    T = cov.filter(pl.col("locus_id").is_in(train_loci)).select(ood_cols).to_numpy().astype(float)
    mu, sd = np.nanmean(T, axis=0), np.nanstd(T, axis=0) + 1e-9
    np.savez(d / "stack" / "ood.npz", mu=mu, sd=sd, T=np.nan_to_num((T - mu) / sd, nan=0.0))
    meta = dict(features=g["features"], transfer_features=t["features"], ood_columns=ood_cols,
                conformal_shift=float(g["conformal_shift"]),
                conformal_shift_by_group={str(k): float(v) for k, v in g["conformal_shift_by_group"].items()},
                group_mode="core", trained_on="212 training people (1000 Genomes), GRCh38 1kGP reference",
                interval_calibration="26 validation people, split-conformal per core-spread group")
    (d / "stack" / "stack.json").write_text(json.dumps(meta, indent=1))

    # ---- score the 52 test people with the package's own code ----
    ms_tmp = d / "loci.parquet"
    L.write_parquet(ms_tmp)
    store.write_manifest(d, name=SET, version="1", source="shipped", reference_fingerprint=fp,
                         reference="GRCh38_1kGP")
    ms = store.ModelSet(d)
    test = sp.filter(pl.col("fold") == "test")["sample_id"].to_list()
    gt = (pl.read_parquet(W / "ground_truth" / "gt_long.parquet", columns=["locus_id", "vamos_sample", "combined_units", "n_haps_called"])
          .filter(pl.col("n_haps_called") == 2)
          .join(pl.DataFrame({"sample_id": list(smap), "vamos_sample": list(smap.values())}), on="vamos_sample"))
    st = ms.stack

    def raw_rows(people):
        out = []
        for sid in people:
            f = pl.read_parquet(W / "features" / f"{sid}.parquet")
            nuis = json.loads((W / "features" / f"{sid}.nuisance.json").read_text())
            for k, v in nuis.items():
                if k not in ("sample_id", "gc_depth_curve"):
                    f = f.with_columns(pl.lit(v, dtype=pl.Float64).alias(f"nuis_{k}"))
            f = add_curve_estimates(f.with_columns(sample_id=pl.lit(sid)), L, single_end=False)
            f = f.join(L.select(["locus_id"] + [c for c in L.columns if c.startswith("cov_")] + ["transfer_ood"]),
                       on="locus_id", how="left")
            out.append(stack.add_transfer(f, st))
        return pl.concat(out, how="diagonal_relaxed")

    # recalibrate the 90% intervals on the validation people (all repeats) with the
    # interval rule that always contains the rounded point estimate
    val = sp.filter(pl.col("fold") == "val")["sample_id"].to_list()
    V = raw_rows(val).join(gt.select("sample_id", "locus_id", "combined_units"), on=["sample_id", "locus_id"])
    lp, lo, hi = stack.raw(V, st)
    t_all, shifts = stack.conformal_shifts(lo, hi, V["combined_units"].to_numpy(), lp, stack.core_groups(V))
    meta["conformal_shift"], meta["conformal_shift_by_group"] = t_all, shifts
    meta["interval_calibration"] = ("26 validation people x all repeats, split-conformal per core-spread group; "
                                    "interval always contains the rounded point estimate")
    (d / "stack" / "stack.json").write_text(json.dumps(meta, indent=1))
    print(f"recalibrated conformal shifts on {V.height:,} validation calls: overall {t_all:+.3f}, {shifts}")
    ms = store.ModelSet(d)
    st = ms.stack
    rows = []
    for sid in test:
        f = stack.predict(raw_rows([sid]), st, t_all, shifts)
        rows.append(tiers(f).select("sample_id", "locus_id", "est_units", "ci90_lo_units", "ci90_hi_units", "tier",
                                    "plc_pred", "transfer_pred"))
    P = pl.concat(rows).join(gt.select("sample_id", "locus_id", "combined_units"), on=["sample_id", "locus_id"])
    perf = report.per_repeat(P)
    L = L.join(perf, on="locus_id", how="left").with_columns(
        perf_source=pl.lit("52 held-out test people (never used in training or calibration)"))
    L.write_parquet(d / "loci.parquet")
    P.write_parquet(args.out / "test_predictions_package.parquet")

    # parity with the research predictions on repeats both scored
    R = pl.read_parquet(W / "model" / "predictions.parquet", columns=["locus_id", "sample_id", "y_pred", "y_lo", "y_hi",
                                                                      "path"]).filter(pl.col("path") == "known")
    J = P.join(R, on=["locus_id", "sample_id"])
    dev = (J["est_units"] - J["y_pred"]).abs()
    print(f"parity vs research on {J.height:,} test calls: point estimates max |diff| {dev.max():.4g}; "
          f"intervals identical {((J['ci90_lo_units'] == J['y_lo']) & (J['ci90_hi_units'] == J['y_hi'])).mean():.4f} "
          f"(differ by design: point-containing intervals, recalibrated shifts)")
    err = (P["est_units"] - P["combined_units"]).abs()
    cover = ((P["combined_units"] >= P["ci90_lo_units"]) & (P["combined_units"] <= P["ci90_hi_units"])).mean()
    print(f"test people, all {L.height:,} repeats: mean abs err {err.mean():.3f}, coverage {cover:.3f}")
    inv = int((P["ci90_lo_units"] > P["ci90_hi_units"]).sum())
    print(f"inverted intervals: {inv}; point estimate outside interval: "
          f"{int(((P['est_units'].round() < P['ci90_lo_units']) | (P['est_units'].round() > P['ci90_hi_units'])).sum())}")
    G = P.join(L.select("locus_id", "cov_panel_core90_units"), on="locus_id")
    G = G.with_columns(group=pl.Series(stack.core_groups(G)), hit=(pl.col("combined_units") >= pl.col("ci90_lo_units"))
                       & (pl.col("combined_units") <= pl.col("ci90_hi_units")))
    print(G.group_by("group").agg(n=pl.len(), coverage=pl.col("hit").mean(),
                                  width=(pl.col("ci90_hi_units") - pl.col("ci90_lo_units")).mean()).sort("group"))

    # ---- TRCompDB annotation lookup for user repeats + default BED ----
    mr = Path("/data/GitHub_Repos/VNTR_Genotyping/data/1000Genomes/mutRate_GRCh38_masterFile.tsv")
    keep = ["chr", "start", "end", "segDupIdentity", "segDupCNV", "consGC", "op1", "op2", "dist2Recom", "size", "period"]
    ann = (pl.read_csv(mr, separator="\t", columns=keep, infer_schema_length=0)
           .rename({"chr": "chrom"}).with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in keep[3:]])
           .with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64)))
    ann.write_parquet(d / "trcompdb_annotation.parquet")
    L.select("chrom", "start", "end", "period", "motifs", "locus_id").sort("chrom", "start").write_csv(
        d / "default_repeats_GRCh38.bed", separator="\t", include_header=False)
    store.write_manifest(d, name=SET, version="1", source="shipped", reference_fingerprint=fp,
                         reference="GRCh38_1kGP (GRCh38_full_analysis_set_plus_decoy_hla)",
                         n_loci=L.height, n_train_people=212, interval_calibration_people=26,
                         performance_people=52, truth="vamos v3.0 (TRCompDB) long-read genotypes",
                         trcompdb_annotation="TRCompDB mutation-rate table v3.0 (Gu & Chaisson)")
    tgz = args.out / "vntrlen-shipped-GRCh38-v1.tar.gz"
    with tarfile.open(tgz, "w:gz") as tf:
        tf.add(d, arcname=f"{fp}/{SET}")
    print(f"wrote {d} and {tgz} ({tgz.stat().st_size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
