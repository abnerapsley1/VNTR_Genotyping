"""
paths.py — repository-root discovery and configuration loading.

Every other module imports `load_config()` and resolves input/output paths
through `Config.path(...)`, so the pipeline can be run from any working
directory.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

_CONFIG_RELPATH = "src/VNTR_AbsoluteGenotyping/config/config.yaml"


def find_repo_root(start: Path | None = None) -> Path:
    """Walk upward from `start` (default: this file) until a `.git` dir is found."""
    here = (start or Path(__file__)).resolve()
    for parent in [here, *here.parents]:
        if (parent / ".git").exists():
            return parent
    raise RuntimeError(
        "Could not locate repository root (no .git found above "
        f"{here}). Pass an explicit config path."
    )


@dataclass
class Config:
    raw: dict[str, Any]
    repo_root: Path
    config_path: Path

    def get(self, dotted: str, default: Any = None) -> Any:
        """Fetch a nested key, e.g. cfg.get('features.mapq_tiers')."""
        node: Any = self.raw
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def path(self, dotted: str, *, must_exist: bool = False) -> Path:
        """Resolve a config value that names a path (relative -> repo root)."""
        val = self.get(dotted)
        if val is None:
            raise KeyError(f"Config key '{dotted}' not found")
        p = Path(val)
        if not p.is_absolute():
            p = self.repo_root / p
        p = p.resolve()
        if must_exist and not p.exists():
            raise FileNotFoundError(f"{dotted} -> {p} does not exist")
        return p

    def ensure_dir(self, dotted: str) -> Path:
        p = self.path(dotted)
        p.mkdir(parents=True, exist_ok=True)
        return p


@lru_cache(maxsize=4)
def load_config(config_path: str | os.PathLike | None = None) -> Config:
    root = find_repo_root()
    cfg_path = Path(config_path).resolve() if config_path else (root / _CONFIG_RELPATH)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")
    with open(cfg_path) as fh:
        raw = yaml.safe_load(fh) or {}
    return Config(raw=raw, repo_root=root, config_path=cfg_path)


if __name__ == "__main__":
    cfg = load_config()
    print(f"repo_root:   {cfg.repo_root}")
    print(f"config_path: {cfg.config_path}")
    print(f"regions_bed: {cfg.path('inputs.regions_bed', must_exist=True)}")
    print(f"vamos_vcf:   {cfg.path('inputs.vamos_vcf', must_exist=True)}")
    print(f"outputs.root:{cfg.path('outputs.root')}")
