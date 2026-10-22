"""Model store: model sets per reference.

    <VNTRLEN_HOME>/models/<reference fingerprint>/<set name>/
        manifest.json    name, version, source (shipped | user), reference, build details
        loci.parquet     one row per repeat: coordinates, motifs, annotation and panel
                         features, calibration curves, performance
        stack/           cross-repeat models (see stack.py), or
        conformal.json   + manifest "stack_from": <set> for user sets that reuse a stack

A repeat's model is identified by (reference fingerprint, locus_id).  When several sets
cover a repeat, the most recently built user set wins unless --model-set is given.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from . import paths
from .stack import Stack


class ModelSet:
    def __init__(self, d: Path):
        self.dir = Path(d)
        self.manifest = json.loads((self.dir / "manifest.json").read_text())
        self._loci = None
        self._stack = None

    name = property(lambda s: s.manifest["name"])
    source = property(lambda s: s.manifest["source"])
    reference_fp = property(lambda s: s.manifest["reference_fingerprint"])

    @property
    def loci(self) -> pl.DataFrame:
        if self._loci is None:
            self._loci = pl.read_parquet(self.dir / "loci.parquet")
        return self._loci

    def stack_dir(self) -> Path:
        src = self.manifest.get("stack_from")
        return self.dir / "stack" if not src else self.dir.parent / src / "stack"

    @property
    def stack(self) -> Stack:
        if self._stack is None:
            self._stack = Stack(self.stack_dir())
        return self._stack

    def conformal(self) -> tuple[float, dict]:
        f = self.dir / "conformal.json"
        src = f if f.exists() else self.stack_dir() / "stack.json"
        c = json.loads(src.read_text())
        return c["conformal_shift"], c["conformal_shift_by_group"]


def all_sets(reference_fp: str | None = None) -> list[ModelSet]:
    out = []
    root = paths.models_dir()
    for m in sorted(root.glob("*/*/manifest.json")):
        ms = ModelSet(m.parent)
        if reference_fp is None or ms.reference_fp == reference_fp:
            out.append(ms)
    return out


def choose(reference_fp: str, locus_ids: list[str], model_set: str | None = None) -> dict[str, ModelSet]:
    """locus_id -> ModelSet that will genotype it (user sets first, newest first)."""
    sets = all_sets(reference_fp)
    if model_set:
        sets = [s for s in sets if s.name == model_set]
        if not sets:
            raise SystemExit(f"model set '{model_set}' not found for this reference")
    sets.sort(key=lambda s: (s.source != "user", -datetime.fromisoformat(s.manifest["created"]).timestamp()))
    want, out = set(locus_ids), {}
    for s in sets:
        for lid in s.loci["locus_id"].to_list():
            if lid in want and lid not in out:
                out[lid] = s
    return out


def new_set_dir(reference_fp: str, name: str) -> Path:
    d = paths.models_dir() / reference_fp / name
    if d.exists():
        raise SystemExit(f"model set '{name}' already exists ({d}); choose another --name")
    d.mkdir(parents=True)
    return d


def write_manifest(d: Path, **kw) -> None:
    kw.setdefault("created", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    (d / "manifest.json").write_text(json.dumps(kw, indent=2, default=str))


QUERY_COLUMNS = ["reference_fingerprint", "model_set", "source", "created", "locus_id", "chrom", "start", "end",
                 "period", "motifs", "n_train_people", "curve_signal", "perf_n_people", "perf_mean_abs_err_units",
                 "perf_within_1_unit", "perf_pearson_r", "perf_ci90_coverage", "perf_source"]


def query(reference_fp: str | None = None, region: str | None = None) -> pl.DataFrame:
    frames = []
    for s in all_sets(reference_fp):
        l = s.loci.with_columns(reference_fingerprint=pl.lit(s.reference_fp), model_set=pl.lit(s.name),
                                source=pl.lit(s.source), created=pl.lit(s.manifest["created"]))
        frames.append(l.select([c for c in QUERY_COLUMNS if c in l.columns]))
    if not frames:
        return pl.DataFrame(schema={c: pl.Utf8 for c in QUERY_COLUMNS})
    df = pl.concat(frames, how="diagonal_relaxed")
    if region:
        chrom, _, rng = region.partition(":")
        df = df.filter(pl.col("chrom") == chrom)
        if rng:
            a, b = (int(x.replace(",", "")) for x in rng.split("-"))
            df = df.filter((pl.col("end") > a - 1) & (pl.col("start") < b))
    return df
