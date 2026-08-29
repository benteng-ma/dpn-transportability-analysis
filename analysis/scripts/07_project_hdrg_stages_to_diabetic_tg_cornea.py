#!/usr/bin/env python3
"""Project frozen human hDRG stage signatures into diabetic TG and cornea datasets."""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
OCULAR_INPUTS = PHASE / "data" / "processed" / "ocular"
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
DATE = "2026-08-27"
SEED = 20260827
N_NULL = 10_000

SIGNATURE_ORDER = [
    "early_allcell_diabetes_vs_control",
    "late_allcell_DPN_vs_diabetes",
    "late_neuron_DPN_vs_diabetes",
    "severity_neuron_modhigh_vs_low_nageotte",
    "xenium_DPN_vs_control",
]
DISPLAY_SIGNATURE = {
    "early_allcell_diabetes_vs_control": "Early all-cell",
    "late_allcell_DPN_vs_diabetes": "Late all-cell",
    "late_neuron_DPN_vs_diabetes": "Late neuron",
    "severity_neuron_modhigh_vs_low_nageotte": "Neuron severity",
    "xenium_DPN_vs_control": "Xenium spatial",
}
DATASET_ORDER = ["GSE329119", "GSE227165", "GSE180490"]
DISPLAY_DATASET = {
    "GSE329119": "Diabetic TG",
    "GSE227165": "Resting diabetic cornea",
    "GSE180490": "Wounded diabetic cornea",
}

DATASETS = {
    "GSE329119": {
        "species": "mouse",
        "tax_id": "10090",
        "compartment": "TG",
        "context": "resting_diabetes",
        "effect": OCULAR_INPUTS / "GSE329119_diabetes_effects.tsv",
        "expression": OCULAR_INPUTS / "GSE329119_gene_FPKM_clean.tsv.gz",
        "effect_gene": "gene_symbol",
        "expression_gene": "gene_symbol",
        "control": ["C1", "C2", "C3"],
        "disease": ["D1", "D2", "D3"],
        "control_label": "control",
        "disease_label": "diabetic",
    },
    "GSE227165": {
        "species": "rat",
        "tax_id": "10116",
        "compartment": "cornea",
        "context": "resting_diabetes",
        "effect": OCULAR_INPUTS / "GSE227165_diabetes_effects_with_mouse_ortholog.tsv",
        "expression": OCULAR_INPUTS / "GSE227165_rat_cornea_gene_FPKM_clean.tsv.gz",
        "effect_gene": "rat_gene_symbol",
        "expression_gene": "gene_name",
        "control": ["NOR_1", "NOR_2", "NOR_3"],
        "disease": ["STZ_1", "STZ_2", "STZ_3"],
        "control_label": "normal",
        "disease_label": "STZ_diabetic",
    },
    "GSE180490": {
        "species": "mouse",
        "tax_id": "10090",
        "compartment": "cornea",
        "context": "wounded_diabetes",
        "effect": OCULAR_INPUTS / "GSE180490_diabetes_effects.tsv",
        "expression": OCULAR_INPUTS / "GSE180490_gene_FPKM_clean.tsv.gz",
        "effect_gene": "gene_symbol",
        "expression_gene": "gene_symbol",
        "control": ["C-Ctrl-W1", "C-Ctrl-W2", "C-Ctrl-W3"],
        "disease": ["C-DM-W1", "C-DM-W2", "C-DM-W3"],
        "control_label": "control_wounded",
        "disease_label": "diabetic_wounded",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = ranked * len(array) / np.arange(1, len(array) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return pd.Series(result, index=values.index)


def load_gene_info(path: Path) -> pd.DataFrame:
    info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)
    return info[["GeneID", "Symbol", "Synonyms"]].copy()


def build_lookup(info: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    exact: dict[str, str] = {}
    folded_candidates: defaultdict[str, set[str]] = defaultdict(set)
    synonym_candidates: defaultdict[str, set[str]] = defaultdict(set)
    id_to_symbol: dict[str, str] = {}
    for row in info.itertuples(index=False):
        gene_id = str(row.GeneID)
        symbol = str(row.Symbol)
        exact[symbol] = gene_id
        folded_candidates[symbol.upper()].add(gene_id)
        id_to_symbol[gene_id] = symbol
        if row.Synonyms and row.Synonyms != "-":
            for synonym in str(row.Synonyms).split("|"):
                if synonym and synonym != "-":
                    synonym_candidates[synonym.upper()].add(gene_id)
    folded = {key: next(iter(ids)) for key, ids in folded_candidates.items() if len(ids) == 1}
    synonyms = {
        key: next(iter(ids))
        for key, ids in synonym_candidates.items()
        if len(ids) == 1 and key not in folded
    }
    return exact, folded, synonyms, id_to_symbol


def resolve(symbol: object, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str]:
    text = str(symbol).strip()
    exact, folded, synonyms, _ = lookup
    if text in exact:
        return exact[text], "official_exact"
    key = text.upper()
    if key in folded:
        return folded[key], "official_casefold"
    if key in synonyms:
        return synonyms[key], "unique_synonym"
    return None, "unresolved_or_ambiguous"


def reciprocal_one_to_one_orthologues(path: Path, other_tax_id: str) -> tuple[dict[str, str], dict[str, int]]:
    pairs: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").split("\t")
        expected = ["#tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"]
        if header != expected:
            raise RuntimeError(f"Unexpected orthologue header: {header}")
        for line in handle:
            tax_id, gene_id, relationship, other_tax, other_gene = line.rstrip("\n").split("\t")
            if relationship != "Ortholog":
                continue
            if tax_id == "9606" and other_tax == other_tax_id:
                pairs.add((gene_id, other_gene))
            elif tax_id == other_tax_id and other_tax == "9606":
                pairs.add((other_gene, gene_id))
    human_to_other: defaultdict[str, set[str]] = defaultdict(set)
    other_to_human: defaultdict[str, set[str]] = defaultdict(set)
    for human, other in pairs:
        human_to_other[human].add(other)
        other_to_human[other].add(human)
    mapping = {
        human: next(iter(others))
        for human, others in human_to_other.items()
        if len(others) == 1 and len(other_to_human[next(iter(others))]) == 1
    }
    return mapping, {
        "unique_pairs": len(pairs),
        "human_with_any": len(human_to_other),
        "other_with_any": len(other_to_human),
        "reciprocal_one_to_one": len(mapping),
    }


def prepare_target(
    dataset: str,
    config: dict[str, object],
    target_lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    effect = pd.read_csv(config["effect"], sep="\t")
    effect["source_row"] = np.arange(2, len(effect) + 2)
    resolved = effect[config["effect_gene"]].map(lambda value: resolve(value, target_lookup))
    effect["target_gene_id"] = [item[0] for item in resolved]
    effect["mapping_method"] = [item[1] for item in resolved]
    _, _, _, id_to_symbol = target_lookup
    effect["target_current_symbol"] = effect["target_gene_id"].map(id_to_symbol)
    for column in ["log2FC", "BH_FDR"]:
        effect[column] = pd.to_numeric(effect[column], errors="coerce")
    if "tested" in effect:
        effect["tested_bool"] = effect["tested"].astype(str).str.lower().isin(["true", "1"])
    else:
        effect["tested_bool"] = False
    effect["fdr_sort"] = effect["BH_FDR"].fillna(np.inf)
    effect["abs_log2fc"] = effect["log2FC"].abs().fillna(-np.inf)
    effect_resolved = effect[effect["target_gene_id"].notna()].sort_values(
        ["target_gene_id", "tested_bool", "fdr_sort", "abs_log2fc", "source_row"],
        ascending=[True, False, True, False, True],
    ).drop_duplicates("target_gene_id", keep="first")

    expression = pd.read_csv(config["expression"], sep="\t", compression="gzip")
    expression["source_row"] = np.arange(2, len(expression) + 2)
    resolved_expr = expression[config["expression_gene"]].map(lambda value: resolve(value, target_lookup))
    expression["target_gene_id"] = [item[0] for item in resolved_expr]
    samples = list(config["control"]) + list(config["disease"])
    for sample in samples:
        expression[sample] = pd.to_numeric(expression[sample], errors="coerce")
    expression["mean_expression"] = expression[samples].mean(axis=1)
    expression["nonmissing_sample_count"] = expression[samples].notna().sum(axis=1)
    expression_resolved = expression[expression["target_gene_id"].notna()].sort_values(
        ["target_gene_id", "nonmissing_sample_count", "mean_expression", "source_row"],
        ascending=[True, False, False, True],
    ).drop_duplicates("target_gene_id", keep="first")
    expression_matrix = expression_resolved.set_index("target_gene_id")[samples]
    expression_mean = expression_resolved.set_index("target_gene_id")["mean_expression"]
    universe = effect_resolved[
        effect_resolved["log2FC"].notna() & effect_resolved["target_gene_id"].isin(expression_matrix.index)
    ].copy()
    universe["mean_expression"] = universe["target_gene_id"].map(expression_mean)
    universe = universe[np.isfinite(universe["mean_expression"]) & (universe["mean_expression"] > 0)].copy()
    expression_matrix = expression_matrix.loc[universe["target_gene_id"]]
    expression_matrix.index.name = "target_gene_id"
    sample_meta = pd.DataFrame(
        {
            "dataset": dataset,
            "sample": samples,
            "group": ["control"] * len(config["control"]) + ["disease"] * len(config["disease"]),
            "source_group": [config["control_label"]] * len(config["control"]) + [config["disease_label"]] * len(config["disease"]),
            "compartment": config["compartment"],
            "context": config["context"],
            "species": config["species"],
        }
    )
    audit = {
        "dataset": dataset,
        "effect_rows": int(len(effect)),
        "resolved_effect_rows": int(effect["target_gene_id"].notna().sum()),
        "unique_resolved_effect_gene_ids": int(len(effect_resolved)),
        "target_universe_gene_ids": int(len(universe)),
        "sample_count": len(samples),
        "effect_sha256": sha256(config["effect"]),
        "expression_sha256": sha256(config["expression"]),
    }
    return universe, expression_matrix, sample_meta, audit


def expression_matched_null(
    up: pd.DataFrame,
    down: pd.DataFrame,
    universe: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    effects = {
        int(decile): frame["log2FC"].to_numpy(dtype=float)
        for decile, frame in universe.groupby("expression_decile")
    }
    up_counts = up["expression_decile"].value_counts().to_dict()
    down_counts = down["expression_decile"].value_counts().to_dict()
    null = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        up_sum = sum(
            float(rng.choice(effects[int(decile)], size=int(count), replace=False).sum())
            for decile, count in up_counts.items()
        )
        down_sum = sum(
            float(rng.choice(effects[int(decile)], size=int(count), replace=False).sum())
            for decile, count in down_counts.items()
        )
        null[iteration] = up_sum / len(up) - down_sum / len(down)
    observed = float(up["log2FC"].mean() - down["log2FC"].mean())
    p = (1 + int(np.sum(null >= observed - 1e-12))) / (N_NULL + 1)
    return null, float(p)


def exact_label_permutation(scores: pd.Series, groups: pd.Series) -> tuple[float, int]:
    n_disease = int((groups == "disease").sum())
    observed = float(scores[groups == "disease"].mean() - scores[groups == "control"].mean())
    values = scores.to_numpy(dtype=float)
    null = []
    for disease_indices in itertools.combinations(range(len(values)), n_disease):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(disease_indices)] = True
        null.append(float(values[mask].mean() - values[~mask].mean()))
    return float(np.mean(np.asarray(null) >= observed - 1e-12)), len(null)


def hedges_g(scores: pd.Series, groups: pd.Series) -> float:
    x = scores[groups == "disease"].to_numpy(dtype=float)
    y = scores[groups == "control"].to_numpy(dtype=float)
    df = len(x) + len(y) - 2
    pooled = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    if pooled <= 0:
        return float("nan")
    correction = 1 - 3 / (4 * df - 1)
    return float(((x.mean() - y.mean()) / math.sqrt(pooled)) * correction)


def leave_one_out(scores: pd.Series, groups: pd.Series) -> tuple[bool, float, float]:
    differences = []
    for sample in scores.index:
        kept = scores.index != sample
        current_scores = scores[kept]
        current_groups = groups[kept]
        differences.append(float(current_scores[current_groups == "disease"].mean() - current_scores[current_groups == "control"].mean()))
    return bool(all(value > 0 for value in differences)), min(differences), max(differences)


def plot_matrix(tests: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    index = tests.set_index(["dataset", "contrast_id"])
    g_matrix = np.array([[index.loc[(dataset, signature), "hedges_g"] for signature in SIGNATURE_ORDER] for dataset in DATASET_ORDER])
    z_matrix = np.array([[
        (index.loc[(dataset, signature), "gene_concordance"] - index.loc[(dataset, signature), "matched_null_mean"])
        / index.loc[(dataset, signature), "matched_null_sd"]
        for signature in SIGNATURE_ORDER
    ] for dataset in DATASET_ORDER])
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 4.8))
    panels = [(g_matrix, "Library-level effect (Hedges g)", -2.5, 2.5), (z_matrix, "Gene concordance vs matched null (SD)", -5, 5)]
    for axis, (matrix, title, vmin, vmax) in zip(axes, panels):
        image = axis.imshow(matrix, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
        axis.set_xticks(range(len(SIGNATURE_ORDER)), [DISPLAY_SIGNATURE[item] for item in SIGNATURE_ORDER], rotation=35, ha="right")
        axis.set_yticks(range(len(DATASET_ORDER)), [DISPLAY_DATASET[item] for item in DATASET_ORDER])
        axis.set_title(title, fontweight="bold")
        for row, dataset in enumerate(DATASET_ORDER):
            for column, signature in enumerate(SIGNATURE_ORDER):
                record = index.loc[(dataset, signature)]
                label = f"{matrix[row, column]:.2f}"
                if bool(record["projection_supportive"]):
                    label += "*"
                axis.text(column, row, label, ha="center", va="center", fontsize=9, color="black")
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    fig.suptitle("Human hDRG disease-stage projection into diabetic trigeminal ganglion and cornea", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.01, "* Frozen provisionally supportive projection; each public dataset has 3 processed libraries per group.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"ocular_hDRG_stage_projection_matrix_{DATE}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    orthologue_path = NCBI / "gene_orthologs.gz"
    info_paths = {
        "human": NCBI / "Homo_sapiens.gene_info.gz",
        "mouse": NCBI / "Mus_musculus.gene_info.complete.gz",
        "rat": NCBI / "Rattus_norvegicus.gene_info.gz",
    }
    lookups = {species: build_lookup(load_gene_info(path)) for species, path in info_paths.items()}
    orthologues = {}
    orthology_qc = {}
    for species, tax_id in [("mouse", "10090"), ("rat", "10116")]:
        orthologues[species], orthology_qc[species] = reciprocal_one_to_one_orthologues(orthologue_path, tax_id)

    signatures = pd.read_csv(SIGNATURES, sep="\t")
    signatures = signatures[signatures["contrast_id"].isin(SIGNATURE_ORDER)].copy()
    resolved = signatures["gene"].map(lambda value: resolve(value, lookups["human"]))
    signatures["human_gene_id"] = [item[0] for item in resolved]
    signatures["human_mapping_method"] = [item[1] for item in resolved]

    all_tests = []
    all_scores = []
    all_mappings = []
    all_sample_meta = []
    dataset_audits = []
    null_summaries = {}
    rng = np.random.default_rng(SEED)
    for dataset in DATASET_ORDER:
        config = DATASETS[dataset]
        universe, expression, sample_meta, audit = prepare_target(dataset, config, lookups[config["species"]])
        universe["expression_decile"] = pd.qcut(universe["mean_expression"].rank(method="first"), 10, labels=False).astype(int)
        universe_index = universe.set_index("target_gene_id")
        target_orthologues = orthologues[config["species"]]
        mapping = signatures.copy()
        mapping["dataset"] = dataset
        mapping["target_species"] = config["species"]
        mapping["target_gene_id"] = mapping["human_gene_id"].map(target_orthologues)
        mapping = mapping.merge(
            universe[["target_gene_id", "target_current_symbol", "log2FC", "BH_FDR", "mean_expression", "expression_decile"]],
            on="target_gene_id",
            how="left",
        )
        mapping["mapped_to_target_universe"] = mapping["log2FC"].notna()
        all_mappings.append(mapping)
        unique_mapping = mapping[mapping["mapped_to_target_universe"]].sort_values(
            ["contrast_id", "direction", "p_val_adj", "gene"]
        ).drop_duplicates(["contrast_id", "direction", "target_gene_id"], keep="first")
        ranks = expression.rank(axis=0, pct=True, method="average") - 0.5
        groups = sample_meta.set_index("sample")["group"].loc[ranks.columns]
        for signature in SIGNATURE_ORDER:
            current = unique_mapping[unique_mapping["contrast_id"] == signature]
            up_ids = current.loc[current["direction"] == "up", "target_gene_id"].tolist()
            down_ids = current.loc[current["direction"] == "down", "target_gene_id"].tolist()
            minimum = 5 if signature == "xenium_DPN_vs_control" else 10
            if len(up_ids) < minimum or len(down_ids) < minimum:
                raise RuntimeError(f"Insufficient mapping for {dataset} / {signature}: up={len(up_ids)}, down={len(down_ids)}")
            up = universe_index.loc[up_ids].reset_index()
            down = universe_index.loc[down_ids].reset_index()
            null, gene_p = expression_matched_null(up, down, universe, rng)
            concordance = float(up["log2FC"].mean() - down["log2FC"].mean())
            scores = ranks.loc[up_ids].mean(axis=0) - ranks.loc[down_ids].mean(axis=0)
            difference = float(scores[groups == "disease"].mean() - scores[groups == "control"].mean())
            exact_p, permutation_count = exact_label_permutation(scores, groups)
            loo_positive, loo_min, loo_max = leave_one_out(scores, groups)
            all_tests.append(
                {
                    "dataset": dataset,
                    "compartment": config["compartment"],
                    "context": config["context"],
                    "species": config["species"],
                    "contrast_id": signature,
                    "mapped_up_n": len(up_ids),
                    "mapped_down_n": len(down_ids),
                    "gene_concordance": concordance,
                    "up_mean_target_log2fc": float(up["log2FC"].mean()),
                    "down_mean_target_log2fc": float(down["log2FC"].mean()),
                    "up_fraction_concordant": float((up["log2FC"] > 0).mean()),
                    "down_fraction_concordant": float((down["log2FC"] < 0).mean()),
                    "matched_null_mean": float(null.mean()),
                    "matched_null_sd": float(null.std(ddof=1)),
                    "gene_matched_one_sided_p": gene_p,
                    "control_score_mean": float(scores[groups == "control"].mean()),
                    "disease_score_mean": float(scores[groups == "disease"].mean()),
                    "disease_minus_control_score": difference,
                    "hedges_g": hedges_g(scores, groups),
                    "exact_one_sided_p": exact_p,
                    "exact_permutation_count": permutation_count,
                    "loo_all_positive": loo_positive,
                    "loo_min_difference": loo_min,
                    "loo_max_difference": loo_max,
                }
            )
            null_summaries[f"{dataset}__{signature}"] = {
                "gene_null_quantiles": {str(q): float(np.quantile(null, q)) for q in [0, 0.025, 0.5, 0.975, 1]}
            }
            for sample, score in scores.items():
                row = sample_meta.set_index("sample").loc[sample]
                all_scores.append(
                    {
                        "dataset": dataset,
                        "contrast_id": signature,
                        "sample": sample,
                        "group": row["group"],
                        "source_group": row["source_group"],
                        "compartment": row["compartment"],
                        "context": row["context"],
                        "species": row["species"],
                        "score": float(score),
                    }
                )
        audit["orthology_species"] = config["species"]
        dataset_audits.append(audit)
        all_sample_meta.append(sample_meta)

    tests = pd.DataFrame(all_tests)
    tests["gene_matched_bh_q"] = bh_adjust(tests["gene_matched_one_sided_p"])
    tests["exact_bh_q"] = bh_adjust(tests["exact_one_sided_p"])
    tests["projection_supportive"] = (
        (tests["disease_minus_control_score"] > 0)
        & (tests["hedges_g"] >= 0.5)
        & tests["loo_all_positive"]
        & (tests["gene_matched_bh_q"] < 0.10)
    )
    test_index = tests.set_index(["dataset", "contrast_id"])
    resting_allcell = any(
        bool(test_index.loc[(dataset, signature), "projection_supportive"])
        for dataset in ["GSE329119", "GSE227165"]
        for signature in ["early_allcell_diabetes_vs_control", "late_allcell_DPN_vs_diabetes"]
    )
    tg_neuronal_signatures = [
        signature
        for signature in ["late_neuron_DPN_vs_diabetes", "severity_neuron_modhigh_vs_low_nageotte"]
        if bool(test_index.loc[("GSE329119", signature), "projection_supportive"])
    ]
    compartment_boundary = any(
        not all(
            bool(test_index.loc[(cornea, signature), "projection_supportive"])
            for cornea in ["GSE227165", "GSE180490"]
        )
        for signature in tg_neuronal_signatures
    )
    ocular_gate = bool(resting_allcell and tg_neuronal_signatures and compartment_boundary)

    tests.to_csv(TABLES / f"ocular_hDRG_stage_projection_tests_{DATE}.tsv", sep="\t", index=False)
    pd.DataFrame(all_scores).to_csv(TABLES / f"ocular_hDRG_stage_projection_library_scores_{DATE}.tsv", sep="\t", index=False)
    pd.concat(all_mappings, ignore_index=True).to_csv(
        TABLES / f"hDRG_to_ocular_one_to_one_signature_mapping_{DATE}.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.concat(all_sample_meta, ignore_index=True).to_csv(
        PHASE / "metadata" / f"ocular_projection_sample_metadata_{DATE}.tsv", sep="\t", index=False
    )
    pd.DataFrame(dataset_audits).to_csv(TABLES / f"ocular_projection_dataset_audit_{DATE}.tsv", sep="\t", index=False)
    plot_matrix(tests)

    qc = {
        "status": "PASS",
        "frozen_spec": "OCULAR_TRANSFER_FROZEN_SPEC_2026-08-27.md",
        "random_seed": SEED,
        "matched_null_iterations": N_NULL,
        "orthology_qc": orthology_qc,
        "input_hashes": {
            str(path): sha256(path)
            for path in [SIGNATURES, orthologue_path, *info_paths.values()]
        },
        "dataset_audits": dataset_audits,
        "null_summaries": null_summaries,
        "gate_components": {
            "resting_ocular_allcell_support": resting_allcell,
            "tg_supportive_neuronal_signatures": tg_neuronal_signatures,
            "minimum_cornea_specificity_boundary": compartment_boundary,
        },
        "ocular_transfer_gate_pass": ocular_gate,
        "interpretation_boundary": (
            "A pass supports anatomy-aware transfer, not propagation, same-animal coupling, temporal causality, or MG/LG diabetes specificity."
        ),
        "excluded_incomplete_download": str(NCBI / "Mus_musculus.gene_info.gz"),
    }
    with (TABLES / f"ocular_hDRG_stage_projection_qc_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    print(tests[[
        "dataset", "contrast_id", "mapped_up_n", "mapped_down_n", "gene_concordance", "gene_matched_bh_q",
        "disease_minus_control_score", "hedges_g", "exact_one_sided_p", "loo_all_positive", "projection_supportive"
    ]].to_string(index=False))
    print(json.dumps({"ocular_transfer_gate_pass": ocular_gate, **qc["gate_components"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
