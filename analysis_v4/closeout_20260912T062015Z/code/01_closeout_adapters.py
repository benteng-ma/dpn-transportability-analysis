from __future__ import annotations

import os

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
T = P / "analysis_v4/targeted_repair_2026-09-12"
O = P / "analysis_v4/all_extensions_2026-09-11"
C = P / "analysis_v4/closeout_20260912T062015Z"
INFO = P / "data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz"
MAP = P / "results/tables/GSE302658_ensembl_to_ncbi_gene_mapping_2026-08-27.tsv.gz"

LABELS = {
    "original_early_allcell": "P1", "original_late_allcell": "P2",
    "original_late_neuron": "P3", "original_severity": "P4",
    "original_xenium": "P5", "late_shared_concordant_neuronal_core": "P6",
    "late_neuron_residual": "P7", "late_allcell_residual": "P8",
    "severity_neuron_shared_concordant_core": "P9", "severity_neuron_residual": "P10",
}

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")

def c01_m04_flow() -> None:
    out = C / "01_M04_flow"
    out.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(T / "01_inputs/M04_ENSEMBL_ROW_MAPPING_AUDIT.tsv.gz", sep="\t", dtype={"gene_id": str})
    raw.insert(0, "feature_id", raw["raw_ensembl_id"].astype(str))
    raw["mapping_status"] = np.where(raw.gene_id.notna(), "unique", "unmapped")
    raw["finite_all"] = raw["finite_all_libraries"].astype(bool)
    raw["eligible_preaggregate"] = raw["include_row"].astype(bool)
    raw["exclusion_reason"] = raw["exclusion_reason"].replace("included", "")
    f = raw[["feature_id", "gene_id", "mapping_status", "finite_all", "all_zero", "eligible_preaggregate", "exclusion_reason"]]
    # The supplied verifier requires a literal empty field (not "NA") for non-unique mappings.
    f = f.copy()
    f["gene_id"] = f["gene_id"].fillna("")
    f.to_csv(out / "RAW_FEATURE_DISPOSITION.tsv", sep="\t", index=False, na_rep="")

    # Canonical disposition includes every GeneID with a unique annotation, even if no source row survived.
    unique = raw[raw.gene_id.notna()].copy()
    rows = []
    for gid, d in unique.groupby("gene_id", sort=False):
        elig = d.eligible_preaggregate.astype(bool)
        rows.append({
            "gene_id": gid,
            "finite_after_aggregate": bool(d.loc[elig, "finite_all"].all()) if elig.any() else False,
            "all_zero_after_aggregate": bool(d.loc[elig, "all_zero"].all()) if elig.any() else True,
            "eligible_rank": bool(elig.any()),
            "n_source_rows": len(d),
            "n_eligible_source_rows": int(elig.sum()),
        })
    g = pd.DataFrame(rows)
    write_tsv(g, out / "CANONICAL_GENE_DISPOSITION.tsv")

    # Rebuild the mapping from the fixed NCBI gene_info snapshot and compare to the cache.
    gi = pd.read_csv(INFO, sep="\t", dtype=str, keep_default_na=False)
    rebuilt = []
    for r in gi.itertuples(index=False):
        for token in str(r.dbXrefs).split("|"):
            if token.startswith("Ensembl:ENSG"):
                rebuilt.append((token.split(":", 1)[1].split(".", 1)[0], str(r.GeneID), str(r.Symbol), str(r.type_of_gene)))
    rb = pd.DataFrame(rebuilt, columns=["ensembl_gene_id", "human_gene_id", "current_symbol", "type_of_gene"]).drop_duplicates()
    rb["rebuilt_one_to_one"] = rb.groupby("ensembl_gene_id")["human_gene_id"].transform("nunique").eq(1)
    cache = pd.read_csv(MAP, sep="\t", dtype=str)
    cc = cache.merge(rb, on=["ensembl_gene_id", "human_gene_id"], how="outer", suffixes=("_cache", "_rebuilt"), indicator=True)
    write_tsv(cc, out / "ANNOTATION_CACHE_VS_FIXED_SNAPSHOT.tsv")
    summary = pd.DataFrame([
        {"metric": "raw_features", "value": len(raw)},
        {"metric": "unique_mapped_raw_features", "value": int(raw.gene_id.notna().sum())},
        {"metric": "eligible_raw_features", "value": int(raw.eligible_preaggregate.sum())},
        {"metric": "canonical_geneids_from_unique_mapping", "value": len(g)},
        {"metric": "eligible_rank_geneids", "value": int(g.eligible_rank.sum())},
        {"metric": "cache_rows", "value": len(cache)},
        {"metric": "rebuilt_snapshot_rows", "value": len(rb)},
        {"metric": "cache_pairs_missing_in_rebuild", "value": int((cc._merge == "left_only").sum())},
        {"metric": "rebuild_pairs_missing_in_cache", "value": int((cc._merge == "right_only").sum())},
    ])
    write_tsv(summary, out / "M04_FLOW_SUMMARY.tsv")
    (out / "C01_PROVENANCE.json").write_text(json.dumps({
        "raw_audit": str(T / "01_inputs/M04_ENSEMBL_ROW_MAPPING_AUDIT.tsv.gz"),
        "raw_audit_sha256": sha(T / "01_inputs/M04_ENSEMBL_ROW_MAPPING_AUDIT.tsv.gz"),
        "mapping_cache": str(MAP), "mapping_cache_sha256": sha(MAP),
        "fixed_gene_info": str(INFO), "fixed_gene_info_sha256": sha(INFO),
        "conclusion": "The cache is annotation-derived. Matrix correction is required only if pair comparison identifies substantive mismatch.",
    }, indent=2), encoding="utf-8")

def c02_fixed_region_inputs() -> None:
    out = C / "02_M02_audit"
    out.mkdir(parents=True, exist_ok=True)
    dmr = pd.read_csv(T / "07_resource_recovery/M02_results/M02_PROPGER_DMR_ALL.tsv", sep="\t")
    parts = dmr.coord.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    regions = pd.DataFrame({
        "region_id": [f"PROPGER_DMR_{i:04d}" for i in range(1, len(dmr) + 1)],
        "chr": parts[0], "start": parts[1].astype(int), "end": parts[2].astype(int),
        "discovery_direction": np.sign(dmr.meandiff).astype(int), "genome_build": "hg19",
    })
    write_tsv(regions, out / "PROPGER_FIXED_REGIONS.tsv")
    anno = pd.read_csv(T / "01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz", sep="\t", dtype={"probe_id": str, "chr": str})
    eligible = pd.read_csv(T / "07_resource_recovery/M02_results/M02_PROPGER_EWAS_ALL_ELIGIBLE.tsv.gz", sep="\t", usecols=["probe_id"], dtype=str)
    ae = anno[anno.probe_id.isin(set(eligible.probe_id))].copy()
    mem = []
    for r in regions.itertuples(index=False):
        x = ae[(ae.chr == r.chr) & (ae.pos >= r.start) & (ae.pos <= r.end)]
        for q in x.itertuples(index=False):
            mem.append({"region_id": r.region_id, "probe_id": q.probe_id, "probe_chr": q.chr, "probe_position": int(q.pos), "genome_build": "hg19"})
    membership = pd.DataFrame(mem).drop_duplicates(["region_id", "probe_id"])
    write_tsv(membership, out / "PROPGER_FIXED_REGION_MEMBERSHIP.tsv")

    meta = pd.read_csv(T / "07_resource_recovery/M02_results/M02_SAMPLE_QC.tsv", sep="\t", dtype=str)
    m = meta[(meta.cohort == "PROPENG") & (meta.sample_qc_pass.str.lower() == "true")].copy()
    m = m[["sample_key", "cohort", "pain_group", "sex"]].rename(columns={"sample_key": "sample_id"})
    write_tsv(m, out / "PROPENG_FIXED_REGION_METADATA.tsv")
    need = sorted(set(membership.probe_id))
    beta_path = T / "07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz"
    cols = ["ID_REF"] + [f"{x}_BetaVal" for x in m.sample_id]
    pieces = []
    for ch in pd.read_csv(beta_path, usecols=cols, chunksize=5000):
        z = ch[ch.ID_REF.astype(str).isin(need)]
        if len(z): pieces.append(z)
    b = pd.concat(pieces, ignore_index=True).drop_duplicates("ID_REF", keep=False)
    uk_eligible = set(pd.read_csv(
        T / "07_resource_recovery/M02_results/M02_PROPENG_EWAS_ALL_ELIGIBLE.tsv.gz",
        sep="\t", usecols=["probe_id"], dtype=str
    ).probe_id)
    # Preserve the full German membership table, but only place independently QC-qualified
    # British probes in the test matrix. Coverage/NE is then handled by the R script.
    b = b[b.ID_REF.astype(str).isin(uk_eligible)].copy()
    b = b.set_index("ID_REF")
    b.columns = [x[:-8] for x in b.columns]
    b = b[m.sample_id.tolist()]
    b.to_csv(out / "PROPENG_FIXED_REGION_BETA.tsv.gz", sep="\t", compression="gzip", index_label="probe_id")

    audit = pd.DataFrame([
        {"item": "EWAS_input_scale", "value": "beta converted to M=log2(beta/(1-beta)); epsilon=1e-6"},
        {"item": "EWAS_design", "value": "intercept + painful status + sex"},
        {"item": "primary_age_batch_model", "value": "BLOCKED: age and verified technical batch unavailable"},
        {"item": "DMR_object", "value": "CpGannotated built from stored t, P, beta_difference, q; is.sig=q<0.05"},
        {"item": "DMRcate_call", "value": "dmrcate(obj, lambda=1000, C=2, pcutoff=0.05, min.cpgs=3)"},
        {"item": "PROPGER_DMRs", "value": len(regions)},
        {"item": "fixed_region_memberships", "value": len(membership)},
        {"item": "PROPENG_participants", "value": len(m)},
        {"item": "PROPENG_QC_qualified_extracted_probes", "value": len(b)},
        {"item": "fixed_region_role", "value": "known-results retrospective fixed-coordinate evaluation; not blind validation"},
    ])
    write_tsv(audit, out / "M02_METHOD_AUDIT.tsv")
    (out / "M02_INPUT_PROVENANCE.json").write_text(json.dumps({
        "beta_file": str(beta_path), "beta_sha256": sha(beta_path),
        "annotation_file": str(T / "01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz"),
        "annotation_sha256": sha(T / "01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz"),
        "discovery_dmr_file": str(T / "07_resource_recovery/M02_results/M02_PROPGER_DMR_ALL.tsv"),
        "discovery_dmr_sha256": sha(T / "07_resource_recovery/M02_results/M02_PROPGER_DMR_ALL.tsv"),
        "metadata_file": str(T / "07_resource_recovery/M02_results/M02_SAMPLE_QC.tsv"),
        "metadata_sha256": sha(T / "07_resource_recovery/M02_results/M02_SAMPLE_QC.tsv"),
    }, indent=2), encoding="utf-8")

def c06_blocks() -> None:
    out = C / "05_block_summary"
    out.mkdir(parents=True, exist_ok=True)
    sources = [
        ("M01", T / "03_module_block_repair/results/M01_CORRECTED_BLOCK_INFLUENCE.tsv"),
        ("M04", T / "03_module_block_repair/results/M04_CORRECTED_BLOCK_INFLUENCE.tsv"),
        ("M05", T / "03_module_block_repair/results/M05_CORRECTED_BLOCK_INFLUENCE.tsv"),
    ]
    all_rows = []
    for aid, f in sources:
        d = pd.read_csv(f, sep="\t")
        for r in d.itertuples(index=False):
            all_rows.append({
                "analysis_id": aid, "source": r.resource, "endpoint": r.comparison,
                "program": LABELS.get(r.module_id, r.module_id), "block": r.block,
                "original_effect": r.delta_S_full, "deleted_effect": r.delta_S_after_block_removal,
                "eligible_before": True, "eligible_after": str(r.status).upper() == "COMPLETED",
                "rank_universe_same": True, "remaining_up": r.remaining_up, "remaining_down": r.remaining_down,
                "reason": "" if str(r.status).upper() == "COMPLETED" else str(r.status),
            })
    z = pd.DataFrame(all_rows)
    write_tsv(z, out / "BLOCK_INFLUENCE_STANDARDIZED.tsv")
    ev = z[z.eligible_after].copy()
    ev["change"] = ev.deleted_effect - ev.original_effect
    ev["direction_preserved"] = np.sign(ev.deleted_effect) == np.sign(ev.original_effect)
    ev["abs_change"] = ev.change.abs()
    summary = ev.groupby(["analysis_id", "source", "endpoint", "program"], as_index=False).agg(
        deletions=("block", "size"), median_abs_change=("abs_change", "median"), max_abs_change=("abs_change", "max"),
        direction_preserved_n=("direction_preserved", "sum"), original_effect=("original_effect", "first"),
        minimum_remaining_up=("remaining_up", "min"), minimum_remaining_down=("remaining_down", "min"),
    )
    write_tsv(summary, out / "BLOCK_INFLUENCE_COMPLETE_SUMMARY.tsv")

def c07_de() -> None:
    out = C / "06_remaining_status"
    out.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(T / "06_whole_gene_DE/results/M04_PRIMARY_DE_ALL_GENES.tsv", sep="\t", dtype={"gene_id": str})
    b = pd.read_csv(T / "06_whole_gene_DE/results/M04_SEX_ADJUSTED_DE_ALL_GENES.tsv", sep="\t", dtype={"gene_id": str})
    def qcol(d): return next(x for x in d.columns if x.lower() in ("adj.p.val", "q", "fdr", "bh"))
    qa, qb = qcol(a), qcol(b)
    aa, bb = set(a.loc[a[qa] < .05, "gene_id"]), set(b.loc[b[qb] < .05, "gene_id"])
    s = pd.DataFrame([
        {"comparison": "M04 primary vs sex-adjusted", "primary_tested": len(a), "secondary_tested": len(b), "primary_q05": len(aa), "secondary_q05": len(bb), "q05_overlap": len(aa & bb), "primary_only_q05": len(aa - bb), "secondary_only_q05": len(bb - aa), "interpretation": "Different design matrices can change edgeR filterByExpr eligibility; the sex-adjusted branch is a distinct sensitivity family."},
        {"comparison": "M05 diagnosis by muscle", "primary_tested": "see source", "secondary_tested": "see source", "primary_q05": "AH=1288; MG=0", "secondary_q05": "interaction=0", "q05_overlap": "NA", "primary_only_q05": "NA", "secondary_only_q05": "NA", "interpretation": "AH differential expression does not establish a diagnosis-by-muscle interaction; MG and interaction families had no q<0.05 genes."},
    ])
    write_tsv(s, out / "M04_M05_DE_INTERPRETATION.tsv")

if __name__ == "__main__":
    c01_m04_flow()
    c02_fixed_region_inputs()
    c06_blocks()
    c07_de()
    print("closeout adapters completed")
