"""vntrlen build — models for a user-defined set of repeats.

Hybrid design:
  * every repeat gets its own calibration curves (fitted on the user's samples);
  * the cross-repeat correction and interval models are reused from the shipped stack,
    with the 90% intervals recalibrated (conformal) on the user's data;
  * with >= 50 people and >= 1,000 repeats (or --retrain always) the correction models
    are retrained on the user's data instead.

Performance is cross-validated over people (5 folds, families kept together): each
person is scored with curves, panel statistics, correction models and interval
calibration that never saw their fold.

Outputs: a new model set in the store, plus in --out:
  performance_per_repeat.tsv, performance_summary.html, build_info.json
"""
from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import polars as pl

from . import (__version__, align, annotate, bed as bedmod, checks, curves, download, extract, gbm, paths,
               report, stack, store, truth as truthmod)

MIN_PEOPLE = 20
WARN_PEOPLE = 50
RETRAIN_PEOPLE, RETRAIN_LOCI = 50, 1000
K = 5


def _fold(sample: str, family: str | None, k: int = K) -> int:
    key = family or sample
    return int.from_bytes(hashlib.sha1(f"vntrlen:{key}".encode()).digest()[:8], "big") % k


def _curve_columns(feats: pl.DataFrame, truth: pl.DataFrame, folds: dict[str, int], single_end: dict[str, bool]):
    """OOF and full-data curves for every repeat."""
    d = feats.join(truth.select("sample_id", "locus_id", "combined_units"), on=["sample_id", "locus_id"])
    d = d.with_columns(fold=pl.col("sample_id").replace_strict(folds, return_dtype=pl.Int64),
                       se=pl.col("sample_id").replace_strict(single_end, return_dtype=pl.Boolean))
    oof_parts, fits = [], {}
    for (lid,), g in d.group_by(["locus_id"], maintain_order=True):
        sig = {c: g[c].cast(pl.Float64).to_numpy() for c in curves.SIGNAL_SHORTLIST if c in g.columns}
        y = g["combined_units"].to_numpy()
        min_n = max(MIN_PEOPLE * (K - 1) // K, 4)
        pred, outside = curves.oof_predict(sig, y, g["fold"].to_numpy(), min_n, g["se"].to_numpy())
        full = curves.fit_locus(sig, y, MIN_PEOPLE)
        fits[lid] = full
        rsd = (full["pe"] or full["se"] or (None, {"resid_sd": np.nan}))[1]["resid_sd"]
        oof_parts.append(g.select("sample_id", "locus_id").with_columns(
            plc_pred=pl.Series(pred), curve_outside_range=pl.Series(outside), plc_resid_sd=pl.lit(rsd)))
    return pl.concat(oof_parts), fits


def run(reference: str | None, bed_path: str, samples: list, truth_path: str, out: Path, name: str,
        annotation: str | None = None, retrain: str = "auto", workers: int = 8, threads: int = 8) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    fasta = download.resolve_reference(reference)
    fp = checks.reference_fingerprint(fasta)
    bed = bedmod.read_bed(bed_path)
    truth = truthmod.read_truth(truth_path, bed)
    have = set(truth["sample_id"])
    samples = [s for s in samples if s.sample_id in have]
    if len(samples) < MIN_PEOPLE:
        raise SystemExit(f"only {len(samples)} samples have ground truth; at least {MIN_PEOPLE} are needed")
    if len(samples) < WARN_PEOPLE:
        print(f"warning: {len(samples)} samples; curves are more reliable with {WARN_PEOPLE}+")
    n_per = truth.group_by("locus_id").len()
    ok_loci = set(n_per.filter(pl.col("len") >= MIN_PEOPLE)["locus_id"])
    dropped = bed.filter(~pl.col("locus_id").is_in(list(ok_loci)))
    bed = bed.filter(pl.col("locus_id").is_in(list(ok_loci)))
    if bed.height == 0:
        raise SystemExit(f"no repeat has ground truth for >= {MIN_PEOPLE} samples")
    print(f"building models for {bed.height:,} repeats with {len(samples)} samples "
          f"({dropped.height:,} repeats dropped: fewer than {MIN_PEOPLE} people with truth)")

    # ---- features ----
    targets = bedmod.to_targets(bed)
    work = paths.cache_dir() / "features"
    frames, single_end = [], {}
    for smp in samples:
        paired = checks.check_short_reads(smp, fasta)
        aln = smp.alignment or align.align_sample(smp, fasta, paths.cache_dir() / "alignments", threads)
        checks.check_alignment_reference(aln, fasta)
        single_end[smp.sample_id] = not paired
        frames.append(extract.extract_sample(smp.sample_id, aln, str(fasta), targets, paired, work, workers))
    feats = pl.concat(frames, how="diagonal_relaxed").filter(pl.col("extraction_ok") == 1.0)

    # ---- annotation, folds, panel statistics, curves ----
    shipped = [s for s in store.all_sets(fp) if s.source == "shipped"] or \
              [s for s in store.all_sets() if s.source == "shipped"]
    if not shipped:
        raise SystemExit("the shipped models are needed as the base stack: vntrlen models install")
    base = shipped[0]
    if base.reference_fp != fp:
        print(f"note: no shipped models for this reference; reusing the GRCh38 cross-repeat models "
              f"(their features are reference-independent ratios, but check performance)")
    trc = base.dir / "trcompdb_annotation.parquet"
    cov = annotate.annotate(bed, trc if trc.exists() else None, Path(annotation) if annotation else None)
    fam = {s.sample_id: s.family_id for s in samples}
    folds = {s: _fold(s, fam.get(s)) for s in fam}
    oof_curve, fits = _curve_columns(feats, truth, folds, single_end)

    st = base.stack
    cov = cov.with_columns(transfer_ood=pl.Series(st.ood_scores(cov)))
    rows = (feats.join(truth, on=["sample_id", "locus_id"])
            .join(oof_curve, on=["sample_id", "locus_id"], how="left")
            .with_columns(plc_available=pl.col("plc_pred").is_finite().cast(pl.Float64),
                          fold=pl.col("sample_id").replace_strict(folds, return_dtype=pl.Int64)))
    # panel statistics for each fold computed without that fold's people
    parts = []
    for f in range(K):
        held = [s for s, k in folds.items() if k == f]
        ps = truthmod.panel_stats(truth, exclude=set(held))
        parts.append(rows.filter(pl.col("fold") == f).join(cov, on="locus_id", how="left")
                     .join(ps, on="locus_id", how="left"))
    rows = stack.add_transfer(pl.concat(parts, how="diagonal_relaxed"), st)

    n_people, n_loci = len(samples), bed.height
    do_retrain = retrain == "always" or (retrain == "auto" and n_people >= RETRAIN_PEOPLE and n_loci >= RETRAIN_LOCI)
    y = rows["combined_units"].to_numpy()
    fold_arr = rows["fold"].to_numpy()
    lp, lo, hi = np.full(y.size, np.nan), np.full(y.size, np.nan), np.full(y.size, np.nan)
    if do_retrain:
        print(f"retraining the correction models on {rows.height:,} rows ({n_people} people x {n_loci} repeats)")
        a = np.log1p(stack.anchor(rows))
        X = stack._X(rows, st.features)
        r = np.log1p(y) - a
        for f in range(K):
            tr, te = fold_arr != f, fold_arr == f
            p = gbm.predict(gbm.train(X[tr], r[tr]), X[te])
            lo[te] = a[te] + np.minimum(gbm.predict(gbm.train(X[tr], r[tr], quantile=ALPHA2), X[te]), p)
            hi[te] = a[te] + np.maximum(gbm.predict(gbm.train(X[tr], r[tr], quantile=1 - ALPHA2), X[te]), p)
            lp[te] = a[te] + p
    else:
        print("reusing the shipped correction models; recalibrating the 90% intervals on your data")
        lp, lo, hi = stack.raw(rows, st)
    groups = stack.core_groups(rows)
    # cross-conformal: each fold's interval uses shifts calibrated on the other folds
    y_lo, y_hi = np.zeros(y.size), np.zeros(y.size)
    for f in range(K):
        te, cal = fold_arr == f, fold_arr != f
        t_all, sh = stack.conformal_shifts(lo[cal], hi[cal], y[cal], lp[cal], groups[cal])
        t = np.array([sh.get(g, t_all) for g in groups[te]])
        y_lo[te], y_hi[te] = stack.interval(lo[te], hi[te], t, lp[te])
    t_all, shifts = stack.conformal_shifts(lo, hi, y, lp, groups)
    rows = rows.with_columns(est_units=pl.Series(np.clip(np.expm1(lp), 0, None)), ci90_lo_units=pl.Series(y_lo),
                             ci90_hi_units=pl.Series(y_hi))

    # ---- save the model set ----
    d = store.new_set_dir(fp, name)
    if do_retrain:
        (d / "stack").mkdir()
        a = np.log1p(stack.anchor(rows))
        X = stack._X(rows, st.features)
        r = np.log1p(y) - a
        gbm.save(gbm.train(X, r), d / "stack" / "point.json")
        gbm.save(gbm.train(X, r, quantile=ALPHA2), d / "stack" / "q_lo.json")
        gbm.save(gbm.train(X, r, quantile=1 - ALPHA2), d / "stack" / "q_hi.json")
        for fname in ("transfer.json", "ood.npz"):
            shutil.copy(st.dir / fname, d / "stack" / fname)
        meta = dict(st.meta, conformal_shift=t_all, conformal_shift_by_group=shifts)
        (d / "stack" / "stack.json").write_text(json.dumps(meta, indent=1))
    else:
        (d / "conformal.json").write_text(json.dumps(dict(conformal_shift=t_all, conformal_shift_by_group=shifts),
                                                     indent=1))
    perf = report.per_repeat(rows)
    npl = dict(n_per.iter_rows())
    full_panel = truthmod.panel_stats(truth)
    bpu = truth.group_by("locus_id").agg(bp_per_unit=pl.col("combined_bp").sum() / pl.col("combined_units").sum())
    loci = (bed.join(cov, on="locus_id", how="left").join(full_panel, on="locus_id", how="left")
            .join(bpu, on="locus_id", how="left").join(perf, on="locus_id", how="left")
            .with_columns(curves=pl.Series([curves.dumps(fits[l]) for l in bed["locus_id"]]),
                          curve_signal=pl.Series([(fits[l]["pe"] or fits[l]["se"] or [None])[0] for l in bed["locus_id"]]),
                          n_train_people=pl.Series([float(npl[l]) for l in bed["locus_id"]]),
                          perf_source=pl.lit("5-fold cross-validation over people")))
    loci.write_parquet(d / "loci.parquet")
    store.write_manifest(d, name=name, version="1", source="user", reference_fingerprint=fp,
                         reference=str(fasta), vntrlen_version=__version__, n_loci=bed.height,
                         n_train_people=n_people, retrained_correction=do_retrain,
                         stack_from=None if do_retrain else base.name, truth_file=str(truth_path),
                         bed=str(bed_path), single_end_samples=sum(single_end.values()))
    perf_full = loci.select("locus_id", "chrom", "start", "end", "n_train_people", "curve_signal",
                            *[c for c in loci.columns if c.startswith("perf_")])
    perf_full.write_csv(out / "performance_per_repeat.tsv", separator="\t", float_precision=4)
    report.html(perf_full, rows, out / "performance_summary.html", name, do_retrain)
    (out / "build_info.json").write_text(json.dumps(dict(
        model_set=name, model_dir=str(d), reference_fingerprint=fp, n_samples=n_people, n_repeats=bed.height,
        repeats_dropped=dropped["locus_id"].to_list()[:1000], retrained_correction=do_retrain,
        created=datetime.now(timezone.utc).isoformat(timespec="seconds")), indent=2))
    print(f"saved model set '{name}' ({bed.height:,} repeats) to {d}; report in {out}")
    return d


ALPHA2 = stack.ALPHA / 2
