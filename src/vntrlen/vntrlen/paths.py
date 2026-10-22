"""Where vntrlen keeps downloaded and user-built data.

    $VNTRLEN_HOME (default ~/.vntrlen)/
        models/<reference fingerprint>/<model set>/   shipped and user-built model sets
        references/                                   downloaded reference FASTA files
        cache/                                        indexes, alignments, scratch
"""
from __future__ import annotations

import os
from pathlib import Path


def home() -> Path:
    h = Path(os.environ.get("VNTRLEN_HOME", Path.home() / ".vntrlen")).expanduser()
    h.mkdir(parents=True, exist_ok=True)
    return h


def models_dir() -> Path:
    d = home() / "models"
    d.mkdir(parents=True, exist_ok=True)
    return d


def references_dir() -> Path:
    d = home() / "references"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cache_dir() -> Path:
    d = home() / "cache"
    d.mkdir(parents=True, exist_ok=True)
    return d
