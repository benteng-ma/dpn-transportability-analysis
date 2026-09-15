#!/usr/bin/env python3
"""DPN closeout helpers. No network access, model selection, or input modification.

Inputs are explicitly adapted TSVs (see SCHEMAS.md), NOT guessed project columns.
Real biological analyses have not been run by the author of this helper package.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable

import numpy as np
import pandas as pd

PROGRAMS = tuple(f"P{i}" for i in range(1, 11))
MODULES = {"black", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def new_directory(path: str | Path) -> Path:
    p = Path(path).resolve()
    require(not p.exists(), f"Refusing to reuse/overwrite output directory: {p}")
    p.mkdir(parents=True)
    return p


def read_tsv(path: str | Path, columns: Iterable[str]) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    require(bool(header) and len(set(header)) == len(header), "Empty or duplicated input header")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require(not df.columns.duplicated().any(), "Duplicate header names are not allowed")
    missing = set(columns) - set(df.columns)
    require(not missing, f"Missing explicitly required columns: {sorted(missing)}")
    return df


def bool_column(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.lower()
    allowed = {"true": True, "false": False, "1": True, "0": False}
    require(values.isin(allowed).all(), f"Non-boolean values in {series.name}")
    return values.map(allowed).astype(bool)


def numeric_column(series: pd.Series, allow_missing: bool = False) -> pd.Series:
    s = series.astype(object).mask(series.astype(str).isin(["", "NA", "NE"]), np.nan)
    x = pd.to_numeric(s, errors="raise")
    require(not np.isinf(x.to_numpy(dtype=float)).any(), f"Infinite input: {series.name}")
    require(allow_missing or not x.isna().any(), f"Missing numeric input: {series.name}")
    return x


def bh_adjust(p: np.ndarray | list[float], family_size: int | None = None) -> np.ndarray:
    """BH on finite P values; NA stays NA. family_size must be explicit when larger."""
    p = np.asarray(p, dtype=float)
    finite = np.isfinite(p)
    require(not np.isinf(p).any(), "Infinite P is invalid")
    require(np.all((p[finite] >= 0) & (p[finite] <= 1)), "P outside [0,1]")
    n = int(finite.sum())
    m = n if family_size is None else int(family_size)
    require(m >= n, "family_size smaller than number of finite tests")
    out = np.full(p.shape, np.nan)
    if not n:
        return out
    idx = np.flatnonzero(finite)
    order = np.argsort(p[idx], kind="stable")
    ranked = p[idx][order]
    adjusted = np.minimum.accumulate((ranked * m / np.arange(1, n + 1))[::-1])[::-1]
    out[idx[order]] = np.minimum(adjusted, 1.0)
    return out


def signed_rank_score(expression: pd.Series, up: Iterable[str], down: Iterable[str],
                      minimum_per_arm: int = 10) -> dict:
    """Score a prevalidated full expression background. This does not choose its rows."""
    require(expression.index.is_unique, "Expression GeneIDs must be unique")
    vals = pd.to_numeric(expression, errors="raise")
    require(np.isfinite(vals.to_numpy()).all(), "Background has non-finite values")
    up, down = set(up), set(down)
    require(not up & down, "A GeneID cannot be in both directions of one program")
    present_up = sorted(up.intersection(vals.index))
    present_down = sorted(down.intersection(vals.index))
    out = {"n_background": len(vals), "n_up": len(present_up), "n_down": len(present_down)}
    if min(len(present_up), len(present_down)) < minimum_per_arm:
        return {**out, "status": "NOT_EVALUABLE", "up_mean": np.nan,
                "down_mean": np.nan, "score": np.nan}
    # Rank ALL background rows before looking up program members.
    r = vals.rank(method="average", pct=True)
    u, d = float(r.loc[present_up].mean()), float(r.loc[present_down].mean())
    return {**out, "status": "COMPLETED", "up_mean": u, "down_mean": d, "score": u - d}


def exact_signflip(differences: np.ndarray, statistic: str) -> tuple[float, float, int]:
    """Full donor sign enumeration, INCLUDING observed signs; no Monte Carlo +1.

    Validity is conditional on the sign-exchangeability/symmetry assumption. This
    is not a randomized-treatment test for observational spatial ROI labels.
    """
    x = np.asarray(differences, dtype=float)
    require(x.ndim == 1 and np.isfinite(x).all(), "Need finite 1D donor differences")
    require(2 <= len(x) <= 16, "Exact helper supports 2..16 donors")
    require(statistic in {"mean", "studentized"}, "Choose the archived statistic explicitly")
    signs = np.array(list(itertools.product((-1.0, 1.0), repeat=len(x))))
    z = signs * x
    means = z.mean(axis=1)
    if statistic == "mean":
        t = np.abs(means)
        observed = abs(float(x.mean()))
    else:
        se = z.std(axis=1, ddof=1) / math.sqrt(len(x))
        t = np.divide(np.abs(means), se, out=np.full(len(means), np.inf), where=se > 0)
        t[(se == 0) & (means == 0)] = 0.0
        observed_se = float(x.std(ddof=1) / math.sqrt(len(x)))
        observed = (abs(float(x.mean())) / observed_se if observed_se > 0
                    else (0.0 if float(x.mean()) == 0 else math.inf))
    if math.isinf(observed):
        exceed = np.isinf(t)
    else:
        tolerance = 64 * np.finfo(float).eps * max(1.0, observed)
        exceed = t >= observed - tolerance
    return float(x.mean()), float(exceed.mean()), len(t)


def stable_rng(seed: int, key: str) -> np.random.Generator:
    digest = hashlib.sha256(f"{seed}|{key}".encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "little"))


def spatial_summary(df: pd.DataFrame, statistic: str, n_boot: int, seed: int,
                    minimum_donors: int) -> pd.DataFrame:
    required = {"donor_id", "program", "effect", "up_effect", "down_effect", "eligible", "reason"}
    require(required.issubset(df.columns), "Incomplete spatial donor schema")
    require(df.program.isin(PROGRAMS).all(), "Unknown program ID")
    require(not df.duplicated(["donor_id", "program"]).any(), "Repeated donor/program")
    require((df.donor_id.str.len() > 0).all(), "Missing donor key")
    work = df.copy()
    work["eligible"] = bool_column(work["eligible"])
    for c in ["effect", "up_effect", "down_effect"]:
        work[c] = numeric_column(work[c], allow_missing=True)
    valid = work.loc[work.eligible]
    require(np.isfinite(valid[["effect", "up_effect", "down_effect"]]).all().all(), "Eligible row lacks effect")
    require(np.allclose(valid.effect, valid.up_effect - valid.down_effect, atol=1e-10, rtol=0),
            "Spatial arm effects do not reconstruct the full effect")
    require(not work.loc[work.program == "P9", "eligible"].any(), "P9 cannot pass the unchanged RNA arm rule")
    require(n_boot >= 100 and minimum_donors >= 2, "Invalid audit resampling configuration")
    rows = []
    for program in PROGRAMS:
        all_rows = work.loc[work.program == program].sort_values("donor_id")
        d = all_rows.loc[all_rows.eligible]
        row = {"program": program, "n_donors": len(d), "statistic": statistic,
               "role": "independent_numerical_audit_not_new_validation"}
        if len(d) < minimum_donors:
            reason = ";".join(sorted(set(all_rows.reason) - {""})) or "no_eligible_donor_records"
            rows.append({**row, "status": "NOT_EVALUABLE", "reason": reason,
                         "effect": np.nan, "P_two_sided": np.nan})
            continue
        x = d.effect.to_numpy(dtype=float)
        effect, p, n = exact_signflip(x, statistic)
        rng = stable_rng(seed, program)
        boot = x[rng.integers(0, len(x), size=(n_boot, len(x)))].mean(axis=1)
        lo, hi = np.quantile(boot, [0.025, 0.975])
        loo = (x.sum() - x) / (len(x) - 1)
        rows.append({**row, "status": "COMPLETED", "reason": "",
                     "effect": effect, "P_two_sided": p, "exact_assignments": n,
                     "audit_ci_low": lo, "audit_ci_high": hi,
                     "audit_bootstrap_draws": n_boot,
                     "delta_up": d.up_effect.mean(), "delta_down": d.down_effect.mean(),
                     "lodo_min_effect": loo.min(), "lodo_max_effect": loo.max(),
                     "lodo_same_sign_fraction": float(np.mean(np.sign(loo) == np.sign(effect))) if effect != 0 else np.nan})
    out = pd.DataFrame(rows)
    out["BH_q_eligible_program_family"] = bh_adjust(out.P_two_sided.to_numpy())
    out["family_finite_n"] = int(out.P_two_sided.notna().sum())
    return out


def flow_audit(features: pd.DataFrame, genes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fcols = {"feature_id", "gene_id", "mapping_status", "finite_all", "all_zero",
             "eligible_preaggregate", "exclusion_reason"}
    gcols = {"gene_id", "finite_after_aggregate", "all_zero_after_aggregate", "eligible_rank",
             "n_source_rows", "n_eligible_source_rows"}
    require(fcols.issubset(features) and gcols.issubset(genes), "Flow schema missing columns")
    f, g = features.copy(), genes.copy()
    require(f.feature_id.is_unique and (f.feature_id.str.len() > 0).all(), "Feature IDs not unique")
    require(g.gene_id.is_unique and (g.gene_id.str.len() > 0).all(), "Gene IDs not unique")
    require(f.mapping_status.isin(["unique", "ambiguous", "unmapped"]).all(), "Invalid mapping status")
    for c in ["finite_all", "all_zero", "eligible_preaggregate"]:
        f[c] = bool_column(f[c])
    for c in ["finite_after_aggregate", "all_zero_after_aggregate", "eligible_rank"]:
        g[c] = bool_column(g[c])
    for c in ["n_source_rows", "n_eligible_source_rows"]:
        x = numeric_column(g[c])
        require(((x >= 0) & (x == np.floor(x))).all(), "Counts must be nonnegative integers")
        g[c] = x.astype(int)
    unique = f.mapping_status == "unique"
    require((f.loc[unique, "gene_id"].str.len() > 0).all(), "Unique mapping lacks GeneID")
    require((f.loc[~unique, "gene_id"] == "").all(), "Ambiguous/unmapped row has guessed GeneID")
    expected_pre = unique & f.finite_all & ~f.all_zero
    require((f.eligible_preaggregate == expected_pre).all(), "Preaggregation mask violates declared strict rule")
    require((f.loc[f.eligible_preaggregate, "exclusion_reason"] == "").all(), "Included feature has exclusion reason")
    require((f.loc[~f.eligible_preaggregate, "exclusion_reason"].str.len() > 0).all(), "Excluded feature lacks reason")
    require(set(f.loc[unique, "gene_id"]) == set(g.gene_id), "Gene table must cover every uniquely mapped gene, including excluded")
    all_counts = f.loc[unique].groupby("gene_id").size()
    pre_counts = f.loc[f.eligible_preaggregate].groupby("gene_id").size()
    require(np.array_equal(g.n_source_rows, g.gene_id.map(all_counts).fillna(0)), "n_source_rows mismatch")
    require(np.array_equal(g.n_eligible_source_rows, g.gene_id.map(pre_counts).fillna(0)), "n_eligible_source_rows mismatch")
    expected_rank = (g.n_eligible_source_rows > 0) & g.finite_after_aggregate & ~g.all_zero_after_aggregate
    require((g.eligible_rank == expected_rank).all(), "Postaggregation eligibility inconsistent")
    n_raw = len(f)
    summary = [
        ("raw_feature", "all", n_raw),
        ("raw_feature", "unique_mapping", int(unique.sum())),
        ("raw_feature", "ambiguous_or_unmapped", int((~unique).sum())),
        ("raw_feature", "all_zero_before_mapping_filter", int(f.all_zero.sum())),
        ("raw_feature", "eligible_preaggregate", int(f.eligible_preaggregate.sum())),
        ("raw_feature", "excluded_preaggregate", int((~f.eligible_preaggregate).sum())),
        ("canonical_gene", "all_uniquely_mapped", len(g)),
        ("canonical_gene", "no_eligible_raw_feature", int((g.n_eligible_source_rows == 0).sum())),
        ("canonical_gene", "eligible_rank", int(g.eligible_rank.sum())),
        ("canonical_gene", "excluded_rank", int((~g.eligible_rank).sum())),
    ]
    reasons = f.loc[~f.eligible_preaggregate].groupby("exclusion_reason").size().rename("n_raw_features").reset_index()
    return pd.DataFrame(summary, columns=["unit", "stage", "count"]), reasons


def block_summary(df: pd.DataFrame, zero_tolerance: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["analysis_id", "source", "endpoint", "program"]
    required = set(keys + ["block", "original_effect", "deleted_effect", "eligible_before",
                          "eligible_after", "rank_universe_same", "remaining_up", "remaining_down", "reason"])
    require(required.issubset(df), "Block schema missing columns")
    d = df.copy()
    require(not d.duplicated(keys + ["block"]).any(), "Duplicated block deletion key")
    require(d.block.isin(MODULES).all(), "Only nine fixed named non-grey modules are valid blocks")
    require(d.program.isin(PROGRAMS).all(), "Unknown program")
    for c in ["eligible_before", "eligible_after", "rank_universe_same"]:
        d[c] = bool_column(d[c])
    for c in ["original_effect", "deleted_effect"]:
        d[c] = numeric_column(d[c], allow_missing=True)
    for c in ["remaining_up", "remaining_down"]:
        d[c] = numeric_column(d[c])
        require(((d[c] >= 0) & (d[c] == np.floor(d[c]))).all(), "Arm counts must be nonnegative integers")
    require(d.rank_universe_same.all(), "Deletion changed rank universe")
    require((~d.eligible_after | d.eligible_before).all(), "Deletion cannot make an originally ineligible program eligible")
    require((d.loc[d.eligible_after, ["remaining_up", "remaining_down"]] >= 10).all().all(), "Deletion arm below 10")
    require(np.isfinite(d.loc[d.eligible_before, "original_effect"]).all(), "Original eligible effect missing")
    require(np.isfinite(d.loc[d.eligible_after, "deleted_effect"]).all(), "Eligible deletion effect missing")
    require((d.loc[~d.eligible_after, "reason"].str.len() > 0).all(), "NE deletion lacks reason")
    d["change"] = np.where(d.eligible_after, d.deleted_effect - d.original_effect, np.nan)
    d["absolute_change"] = d.change.abs()
    d["baseline_near_zero"] = d.original_effect.abs() <= zero_tolerance
    d["signed_retention_ratio"] = np.nan
    ix = d.eligible_after & ~d.baseline_near_zero
    d.loc[ix, "signed_retention_ratio"] = d.loc[ix, "deleted_effect"] / d.loc[ix, "original_effect"]
    d["same_sign"] = np.nan
    d.loc[ix, "same_sign"] = (np.sign(d.loc[ix, "deleted_effect"]) == np.sign(d.loc[ix, "original_effect"])).astype(float)
    rows = []
    for key, group in d.groupby(keys, sort=True, dropna=False):
        before = group.loc[group.eligible_before, "original_effect"]
        if len(before):
            require(np.allclose(before, before.iloc[0], rtol=0, atol=1e-10), "Baseline effect differs across deletions")
        finite = group.loc[group.eligible_after]
        row = dict(zip(keys, key))
        row.update(n_planned=len(group), n_evaluable=len(finite), n_ne=int((~group.eligible_after).sum()))
        if len(finite):
            m = float(finite.absolute_change.max())
            row.update(max_absolute_change=m, median_absolute_change=float(finite.absolute_change.median()),
                       same_sign_fraction=float(finite.same_sign.mean()) if finite.same_sign.notna().any() else np.nan,
                       max_change_blocks=";".join(sorted(finite.loc[np.isclose(finite.absolute_change, m, rtol=0, atol=1e-12), "block"])),
                       status="COMPLETED_DESCRIPTIVE")
        else:
            row.update(status="NOT_EVALUABLE", max_absolute_change=np.nan, median_absolute_change=np.nan,
                       same_sign_fraction=np.nan, max_change_blocks="")
        rows.append(row)
    return pd.DataFrame(rows), d


def inventory(root: Path, output: Path) -> None:
    require(root.name == "phase0_6_human_dpn_stage_projection", "Point --root at the actual DPN project root")
    skip = {".git", ".venv", "venv", "node_modules", "renv", "__pycache__", "site-library", "packrat"}
    patterns = ("locked", "status", "dmr", "roi", "block", "background", "manifest", "membership", "scale", "sessioninfo")
    rows = []
    for f in sorted(root.rglob("*")):
        rel = f.relative_to(root)
        if any(x in skip for x in rel.parts) or f.is_symlink() or not f.is_file():
            continue
        if f.resolve().is_relative_to(output):
            continue
        if any(p in f.name.lower() for p in patterns):
            rows.append({"path": rel.as_posix(), "bytes": f.stat().st_size,
                         "role": "CANDIDATE_ONLY_NOT_SELECTED", "mtime_ns": f.stat().st_mtime_ns})
    pd.DataFrame(rows, columns=["path", "bytes", "role", "mtime_ns"]).to_csv(output / "CANDIDATE_FILES.tsv", sep="\t", index=False)


def tar_inventory(path: Path, output: Path) -> None:
    rows = []
    with tarfile.open(path, "r:*") as archive:
        for entry in archive:
            pure = PurePosixPath(entry.name)
            safe = not pure.is_absolute() and ".." not in pure.parts and not entry.issym() and not entry.islnk()
            name = entry.name.lower()
            role = ("SCALE_CANDIDATE" if "scalefactor" in name else
                    "IMAGE_CANDIDATE" if name.endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff")) else
                    "ROI_CANDIDATE" if "roi" in name or "selection" in name else "OTHER")
            rows.append({"member": entry.name, "bytes": entry.size, "is_regular_file": entry.isfile(),
                         "safe_path": safe, "candidate_role": role})
    pd.DataFrame(rows).to_csv(output / "TAR_MEMBERS.tsv", sep="\t", index=False)
    require(all(r["safe_path"] for r in rows), "Unsafe archive path/link detected; no extraction performed")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    for name in ["flow", "blocks", "spatial", "inventory", "tar-inventory"]:
        p = sub.add_parser(name)
        p.add_argument("--out", required=True, help="A NEW directory, not a historical result directory")
        if name == "flow":
            p.add_argument("--features", required=True)
            p.add_argument("--genes", required=True)
        elif name == "blocks":
            p.add_argument("--input", required=True)
            p.add_argument("--zero-tolerance", type=float, required=True)
        elif name == "spatial":
            p.add_argument("--input", required=True)
            p.add_argument("--statistic", choices=["mean", "studentized"], required=True)
            p.add_argument("--bootstrap-draws", type=int, required=True)
            p.add_argument("--seed", type=int, required=True)
            p.add_argument("--minimum-donors", type=int, required=True)
        elif name == "inventory":
            p.add_argument("--root", required=True)
        else:
            p.add_argument("--tar", required=True)
    a = ap.parse_args()
    output = new_directory(a.out)
    start = datetime.now(timezone.utc).isoformat()
    receipt = {"command": a.command, "started_utc": start, "status": "RUNNING",
               "claim": "helper execution; completion does not prove clinical validation"}
    try:
        if a.command == "flow":
            s, r = flow_audit(read_tsv(a.features, []), read_tsv(a.genes, []))
            s.to_csv(output / "COUNT_FLOW.tsv", sep="\t", index=False)
            r.to_csv(output / "RAW_EXCLUSIONS.tsv", sep="\t", index=False)
        elif a.command == "blocks":
            require(a.zero_tolerance >= 0, "Negative tolerance")
            s, r = block_summary(read_tsv(a.input, []), a.zero_tolerance)
            s.to_csv(output / "BLOCK_SUMMARY.tsv", sep="\t", index=False, na_rep="NA")
            r.to_csv(output / "BLOCK_DETAIL.tsv", sep="\t", index=False, na_rep="NA")
        elif a.command == "spatial":
            s = spatial_summary(read_tsv(a.input, []), a.statistic, a.bootstrap_draws, a.seed, a.minimum_donors)
            s.to_csv(output / "SPATIAL_AUDIT.tsv", sep="\t", index=False, na_rep="NA")
        elif a.command == "inventory":
            inventory(Path(a.root).resolve(), output)
        else:
            tar_inventory(Path(a.tar).resolve(), output)
        receipt["status"] = "COMPLETED"
    except Exception as exc:
        receipt.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        receipt["ended_utc"] = datetime.now(timezone.utc).isoformat()
        (output / "HELPER_RECEIPT.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
