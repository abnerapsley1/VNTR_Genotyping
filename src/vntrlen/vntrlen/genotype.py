"""vntrlen genotype — combined repeat length (both alleles) per sample and repeat.

Outputs in --out:
  genotypes.tsv            one row per sample x repeat (separate sample_id and locus_id columns)
  genotypes_matrix.tsv     repeats x samples, point estimates in repeat units
  skipped_repeats.tsv      BED repeats without a model for this reference
  run_info.json            inputs, reference fingerprint, model sets, versions
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl

from . import __version__, align, bed as bedmod, checks, curves, download, extract, paths, stack, store

OUT_COLUMNS = ["sample_id", "locus_id", "chrom", "start", "end", "period", "est_units", "est_bp",
               "ci90_lo_units", "ci90_hi_units", "ci90_lo_bp", "ci90_hi_bp", "tier", "flags",
               "model_set", "model_version"]


def tiers(df: pl.DataFrame) -> pl.DataFrame:
    width_rel = (pl.col("ci90_hi_units") - pl.col("ci90_lo_units")) / pl.max_horizontal(pl.col("est_units"), pl.lit(1.0))
    return df.with_columns(tier=pl.when(~pl.col("est_units").is_finite() | (pl.col("extraction_ok") != 1.0))
                           .then(pl.lit("NA"))
                           .when((width_rel <= 1.0) & ~pl.col("curve_outside_range")).then(pl.lit("PASS"))
                           .otherwise(pl.lit("LOW_CONF")))


_PARSED: dict[int, dict] = {}


def parsed_curves(loci: pl.DataFrame) -> dict:
    """locus_id -> parsed curves, cached per loci table."""
    key = id(loci)
    if key not in _PARSED:
        _PARSED[key] = {l: curves.loads(c) for l, c in zip(loci["locus_id"].to_list(), loci["curves"].to_list())}
    return _PARSED[key]


def add_curve_estimates(feat: pl.DataFrame, loci: pl.DataFrame, single_end: bool) -> pl.DataFrame:
    """plc_pred / plc_resid_sd / plc_available / curve_outside_range from each repeat's stored curves."""
    pc = parsed_curves(loci)
    n = feat.height
    pred, rsd, out = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, dtype=bool)
    cols = {c: feat[c].cast(pl.Float64).to_numpy() for c in curves.SIGNAL_SHORTLIST if c in feat.columns}
    for i, lid in enumerate(feat["locus_id"].to_list()):
        c = pc.get(lid)
        if not c:
            continue
        entry = c.get("se" if single_end else "pe") or c.get("pe") or c.get("se")
        if not entry:
            continue
        rsd[i] = entry[1]["resid_sd"]     # the repeat's curve reliability, even if this person's signal is missing
        pv, ov = curves.predict_entry(entry, {k: v[i:i + 1] for k, v in cols.items()})
        if np.isfinite(pv[0]):
            pred[i], out[i] = pv[0], bool(ov[0])
    return feat.with_columns(plc_pred=pl.Series(pred), plc_resid_sd=pl.Series(rsd),
                             plc_available=pl.Series(np.isfinite(pred).astype(float)),
                             curve_outside_range=pl.Series(out))


def predict_rows(feat: pl.DataFrame, ms: store.ModelSet, single_end: bool) -> pl.DataFrame:
    loci = ms.loci
    cov_cols = [c for c in loci.columns if c.startswith("cov_") or c in ("transfer_ood", "bp_per_unit")]
    df = add_curve_estimates(feat, loci, single_end).join(
        loci.select(["locus_id", "chrom", "start", "end", "period"] + cov_cols), on="locus_id", how="left",
        suffix="_bed")
    st = ms.stack
    df = stack.add_transfer(df, st)
    t_all, shifts = ms.conformal()
    df = stack.predict(df, st, t_all, shifts)
    bpu = pl.col("bp_per_unit")
    df = df.with_columns(est_bp=pl.col("est_units") * bpu, ci90_lo_bp=pl.col("ci90_lo_units") * bpu,
                         ci90_hi_bp=pl.col("ci90_hi_units") * bpu)
    df = tiers(df)
    flags = []
    for oob in df["curve_outside_range"].to_list():
        f = []
        if single_end:
            f.append("single_end_unvalidated")
        if oob:
            f.append("signal_outside_training_range")
        flags.append(",".join(f) or ".")
    return df.with_columns(flags=pl.Series(flags), model_set=pl.lit(ms.name),
                           model_version=pl.lit(ms.manifest.get("version", "?")))


def run(reference: str | None, bed_path: str | None, samples: list, out: Path, model_set: str | None = None,
        workers: int = 8, threads: int = 8) -> pl.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    fasta = download.resolve_reference(reference)
    fp = checks.reference_fingerprint(fasta)
    sets = store.all_sets(fp)
    if not sets:
        raise SystemExit(f"no models installed for this reference (fingerprint {fp}). For the shipped GRCh38 "
                         "models run: vntrlen models install")
    if bed_path:
        bed = bedmod.read_bed(bed_path)
    else:
        shipped = [s for s in sets if s.source == "shipped"]
        if not shipped:
            raise SystemExit("no BED given and no shipped model set for this reference")
        bed = shipped[0].loci.select(list(bedmod.SCHEMA) + ["locus_id"])
    chosen = store.choose(fp, bed["locus_id"].to_list(), model_set)
    skipped = bed.filter(~pl.col("locus_id").is_in(list(chosen)))
    skipped.select("locus_id", "chrom", "start", "end").with_columns(
        reason=pl.lit("no model for this reference; build one with: vntrlen build")).write_csv(
        out / "skipped_repeats.tsv", separator="\t")
    if skipped.height:
        print(f"{skipped.height:,} of {bed.height:,} repeats have no model and are skipped "
              f"(see {out / 'skipped_repeats.tsv'})")
    by_set: dict[str, list[str]] = {}
    for lid, s in chosen.items():
        by_set.setdefault(s.name, []).append(lid)
    set_obj = {s.name: s for s in sets}
    work = paths.cache_dir() / "features"
    results = []
    for smp in samples:
        paired = checks.check_short_reads(smp, fasta)
        aln = smp.alignment
        if smp.needs_alignment:
            aln = align.align_sample(smp, fasta, paths.cache_dir() / "alignments", threads)
        checks.check_alignment_reference(aln, fasta)
        for name, lids in by_set.items():
            ms = set_obj[name]
            # extract with the model set's own motifs, so features match how the model was built
            tg = bedmod.to_targets(ms.loci.filter(pl.col("locus_id").is_in(lids)))
            feat = extract.extract_sample(smp.sample_id, aln, str(fasta), tg, paired, work, workers)
            results.append(predict_rows(feat, ms, single_end=not paired).select(OUT_COLUMNS))
    res = pl.concat(results).sort("sample_id", "chrom", "start").with_columns(
        pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64), pl.col("period").cast(pl.Int64))
    res.write_csv(out / "genotypes.tsv", separator="\t", float_precision=3)
    res.pivot(on="sample_id", index="locus_id", values="est_units").write_csv(
        out / "genotypes_matrix.tsv", separator="\t", float_precision=3)
    (out / "run_info.json").write_text(json.dumps(dict(
        vntrlen_version=__version__, reference=str(fasta), reference_fingerprint=fp,
        bed=bed_path or "shipped default", n_repeats_genotyped=len(chosen), n_repeats_skipped=skipped.height,
        samples=[s.sample_id for s in samples], model_sets=sorted(by_set)), indent=2))
    print(f"wrote {out / 'genotypes.tsv'} ({res.height:,} rows)")
    return res
