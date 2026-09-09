"""
manifest.py — read the user-filled config/samples.tsv and config/references.tsv.

Both are tab-separated; '#' comment lines and blank lines are ignored.  Column
order is fixed by the templates; a header row (first non-comment line whose first
field is 'sample_id' / 'ref_key') is optional and skipped if present.
"""

from __future__ import annotations

from pathlib import Path

_SAMPLE_FIELDS = ["sample_id", "vamos_column", "cohort", "data_url",
                  "index_url", "align_reference", "md5"]
_REF_FIELDS = ["ref_key", "fasta_url", "fai_url", "contig_style", "build", "notes"]


def _read_rows(path: Path, fields: list[str], header_key: str) -> list[dict]:
    if not path.exists():
        tmpl = path.with_name(path.stem + ".template" + path.suffix)
        hint = f"  (copy {tmpl.name} -> {path.name} and fill it in)" if tmpl.exists() else ""
        raise SystemExit(f"manifest not found: {path}{hint}")
    rows: list[dict] = []
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        parts = ln.split("\t")
        if parts[0].strip() == header_key:
            continue
        parts = (parts + [""] * len(fields))[: len(fields)]
        row = {k: (v.strip() or None) for k, v in zip(fields, parts)}
        rows.append(row)
    return rows


def load_samples(cfg) -> list[dict]:
    rows = _read_rows(cfg.path("inputs.samples_tsv"), _SAMPLE_FIELDS, "sample_id")
    if not rows:
        raise SystemExit(
            f"No sample rows in {cfg.path('inputs.samples_tsv')} — fill the manifest."
        )
    seen = set()
    for r in rows:
        if not r["sample_id"] or not r["data_url"]:
            raise SystemExit(f"manifest row missing sample_id/data_url: {r}")
        if r["sample_id"] in seen:
            raise SystemExit(f"duplicate sample_id in manifest: {r['sample_id']}")
        seen.add(r["sample_id"])
    return rows


def load_references(cfg) -> dict[str, dict]:
    path = cfg.path("inputs.references_tsv")
    if not path.exists():
        return {}
    return {r["ref_key"]: r for r in _read_rows(path, _REF_FIELDS, "ref_key")
            if r["ref_key"]}
