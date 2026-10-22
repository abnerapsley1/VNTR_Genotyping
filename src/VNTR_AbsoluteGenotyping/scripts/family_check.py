"""Family checks of per-haplotype calls on held-out (test) families.

Trio SH032 (father HG00512, mother HG00513, child HG00514) and duo GB40
(father HG02666, son HG02668) are all in the test fold, so none of them
informed the model or the allele catalog.

  Mendelian consistency  child's two alleles = one from each parent (exact
                         whole repeat units), at known loci where all members
                         have a call; compared for vamos truth, the catalog
                         method, the SNP-haplotype method, and equal halves.
  Duo                    child shares at least one allele with the father.
  Phase (SNP method)     the child's 1kGP haplotypes are pedigree-phased, so
                         each SNP haplotype comes from one parent along a whole
                         chromosome.  Where only one parent carries the
                         child's longer allele, that names the haplotype that
                         should carry it; agreement with the SNP method's
                         phase call is reported per chromosome as
                         max(f, 1 - f) (the panel's hap1/hap2 order may be
                         paternal-first or maternal-first).

    python scripts/family_check.py  -> work/model/family_check.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.model.assemble import _clean_vamos_col  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402

TRIO = {"father": "HG00512", "mother": "HG00513", "child": "HG00514"}
DUO = {"father": "HG02666", "child": "HG02668"}


def calls(cfg) -> pl.DataFrame:
    mdir = cfg.path("outputs.model_dir")
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg) if s.get("vamos_column")}
    gt = (pl.read_parquet(cfg.path("outputs.ground_truth_dir") / "gt_long.parquet",
                          columns=["locus_id", "vamos_sample", "a1_units", "a2_units", "n_haps_called"])
          .filter(pl.col("n_haps_called") == 2)
          .join(pl.DataFrame({"sample_id": list(smap), "vamos_sample": list(smap.values())}), on="vamos_sample")
          .select("locus_id", "sample_id", t_lo=pl.min_horizontal("a1_units", "a2_units"),
                  t_hi=pl.max_horizontal("a1_units", "a2_units")))
    h = pl.read_parquet(mdir / "haplotype_predictions.parquet",
                        columns=["locus_id", "sample_id", "path", "h_lo", "h_hi", "half_lo", "half_hi"])
    d = h.filter(pl.col("path") == "known").join(gt, on=["locus_id", "sample_id"], how="left")
    f = mdir / "haplotype_snp_predictions.parquet"
    if f.exists():
        have = pl.read_parquet_schema(f)
        tag = "full" if "full_lo" in have else "snp"           # best available SNP-based variant
        sel = [pl.col("locus_id"), pl.col("sample_id"), pl.col(f"{tag}_lo").alias("s_lo"),
               pl.col(f"{tag}_hi").alias("s_hi"), pl.col(f"{tag}_phase_hap1_long").alias("s_phase_hap1_long")]
        if tag == "full":
            sel += [pl.col("snp_lo"), pl.col("snp_hi"), pl.col("cat_span_lo").alias("cs_lo"),
                    pl.col("cat_span_hi").alias("cs_hi")]
        s = pl.read_parquet(f).select(sel)
        d = d.join(s, on=["locus_id", "sample_id"], how="left")
    return d


def mendel(c, m, f) -> np.ndarray:
    """c, m, f: (n, 2) allele pairs -> child pair explainable by one allele from each parent."""
    def has(p, x):
        return (p[:, 0] == x) | (p[:, 1] == x)
    return (has(m, c[:, 0]) & has(f, c[:, 1])) | (has(m, c[:, 1]) & has(f, c[:, 0]))


def main():
    cfg = load_config(None)
    d = calls(cfg)
    methods = [("truth (vamos)", "t_lo", "t_hi"), ("catalog method", "h_lo", "h_hi"), ("equal halves", "half_lo", "half_hi")]
    if "snp_lo" in d.columns and "s_lo" in d.columns:
        methods.insert(2, ("SNP haplotypes", "snp_lo", "snp_hi"))
    if "cs_lo" in d.columns:
        methods.insert(2, ("catalog + spanning reads", "cs_lo", "cs_hi"))
    if "s_lo" in d.columns:
        methods.insert(2, ("best SNP-based method", "s_lo", "s_hi"))
    rep = {"trio": TRIO, "duo": DUO, "mendelian_consistency": {}, "duo_shares_allele": {}}
    roles = {r: d.filter(pl.col("sample_id") == s).drop("sample_id", "path") for r, s in TRIO.items()}
    t = roles["child"].join(roles["mother"], on="locus_id", suffix="_m").join(roles["father"], on="locus_id", suffix="_f")
    for name, lo, hi in methods:
        cols = [lo, hi, f"{lo}_m", f"{hi}_m", f"{lo}_f", f"{hi}_f"]
        x = t.drop_nulls(cols)
        ok = mendel(x.select(lo, hi).to_numpy(), x.select(f"{lo}_m", f"{hi}_m").to_numpy(),
                    x.select(f"{lo}_f", f"{hi}_f").to_numpy())
        # informative loci: truth shows the child's two alleles differ
        inf = (x["t_lo"] != x["t_hi"]).to_numpy() if "t_lo" in x.columns else np.zeros(x.height, bool)
        rep["mendelian_consistency"][name] = {"n_loci": x.height, "consistent": float(ok.mean()),
                                              "consistent_where_child_het": float(ok[inf].mean()) if inf.any() else None}
        print(f"trio {name:22s}: {ok.mean():.3f} of {x.height:,} loci Mendelian-consistent"
              + (f" ({ok[inf].mean():.3f} where the child's alleles differ)" if inf.any() else ""))
    droles = {r: d.filter(pl.col("sample_id") == s).drop("sample_id", "path") for r, s in DUO.items()}
    du = droles["child"].join(droles["father"], on="locus_id", suffix="_f")
    for name, lo, hi in methods:
        x = du.drop_nulls([lo, hi, f"{lo}_f", f"{hi}_f"])
        c, f = x.select(lo, hi).to_numpy(), x.select(f"{lo}_f", f"{hi}_f").to_numpy()
        ok = (c[:, :, None] == f[:, None, :]).any((1, 2))
        rep["duo_shares_allele"][name] = {"n_loci": x.height, "shares": float(ok.mean())}
        print(f"duo  {name:22s}: {ok.mean():.3f} of {x.height:,} loci share an allele")

    if "s_phase_hap1_long" in t.columns:
        # loci where truth identifies the parent of the child's longer allele (only one parent carries it,
        # and the other parent carries the child's shorter allele)
        x = t.drop_nulls(["t_lo", "t_hi", "t_lo_m", "t_hi_m", "t_lo_f", "t_hi_f", "s_phase_hap1_long", "s_lo", "s_hi"])
        x = x.filter((pl.col("t_lo") != pl.col("t_hi")) & (pl.col("s_lo") == pl.col("t_lo")) & (pl.col("s_hi") == pl.col("t_hi")))
        inm = (pl.col("t_hi") == pl.col("t_lo_m")) | (pl.col("t_hi") == pl.col("t_hi_m"))
        inf_ = (pl.col("t_hi") == pl.col("t_lo_f")) | (pl.col("t_hi") == pl.col("t_hi_f"))
        x = x.with_columns(long_from=pl.when(inm & ~inf_).then(pl.lit("mother")).when(inf_ & ~inm).then(pl.lit("father")))
        x = x.drop_nulls("long_from").with_columns(chrom=pl.col("locus_id").str.split(":").list.get(0))
        x = x.filter((pl.col("s_phase_hap1_long") - 0.5).abs() >= 0.3)        # confident phase calls only
        per = []
        for (ch,), g in x.group_by(["chrom"], maintain_order=True):
            hap1_long = (g["s_phase_hap1_long"] > 0.5).to_numpy()
            from_father = (g["long_from"] == "father").to_numpy()
            agree = float(np.mean(hap1_long == from_father))
            per.append({"chrom": ch, "n": g.height, "agreement": max(agree, 1 - agree),
                        "hap1_is_paternal": agree >= 0.5})
        n = sum(p["n"] for p in per)
        overall = sum(p["agreement"] * p["n"] for p in per) / n if n else None
        rep["child_phase"] = {"n_loci": n, "agreement": overall, "per_chromosome": per}
        print(f"child phase (SNP method, confident calls): {overall:.3f} agreement over {n:,} informative loci"
              if n else "child phase: no informative loci")
    out = cfg.path("outputs.model_dir") / "family_check.json"
    out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
