#!/usr/bin/env python3
"""Build a frozen, explanatory cross-target atlas for source-defined hDRG components."""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
SCRIPTS = PHASE / "analysis" / "scripts"
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
DATE = "2026-08-27"
SEED = 20260827
N_NULL = 10_000
MODULE_IDS = [
    "late_shared_concordant_neuronal_core",
    "late_neuron_residual",
    "late_allcell_residual",
]
MODULE_DISPLAY = {
    "late_shared_concordant_neuronal_core": "Shared neuronal core",
    "late_neuron_residual": "Neuron residual",
    "late_allcell_residual": "All-cell residual",
}


def load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


M12 = load_script("phase06_sural", "12_validate_hdrg_components_in_human_sural_nerve.py")
M06 = load_script("phase06_human_bulk", "06_validate_hdrg_signatures_in_independent_human_bulk.py")
M07 = load_script("phase06_ocular", "07_project_hdrg_stages_to_diabetic_tg_cornea.py")


def exact_permutation(
    scores: pd.Series,
    groups: pd.Series,
    positive: str,
    negative: str,
    two_sided: bool,
) -> tuple[float, int]:
    keep = groups.isin([positive, negative])
    scores = scores[keep]
    groups = groups[keep]
    n_positive = int((groups == positive).sum())
    observed = float(scores[groups == positive].mean() - scores[groups == negative].mean())
    values = scores.to_numpy(dtype=float)
    null = []
    for positive_indices in itertools.combinations(range(len(values)), n_positive):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(positive_indices)] = True
        null.append(float(values[mask].mean() - values[~mask].mean()))
    null = np.asarray(null)
    if two_sided:
        p_value = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
    else:
        p_value = float(np.mean(null >= observed - 1e-12))
    return p_value, len(null)


def hedges_g(scores: pd.Series, groups: pd.Series, positive: str, negative: str) -> float:
    x = scores[groups == positive].to_numpy(dtype=float)
    y = scores[groups == negative].to_numpy(dtype=float)
    df = len(x) + len(y) - 2
    if len(x) < 2 or len(y) < 2 or df <= 0:
        return np.nan
    pooled = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    if pooled <= 0:
        return np.nan
    correction = 1 - 3 / (4 * df - 1)
    return float(((x.mean() - y.mean()) / math.sqrt(pooled)) * correction)


def leave_one_out(
    scores: pd.Series, groups: pd.Series, positive: str, negative: str
) -> tuple[bool, float, float]:
    keep = groups.isin([positive, negative])
    scores = scores[keep]
    groups = groups[keep]
    differences = []
    for sample in scores.index:
        retained = scores.index != sample
        current_scores = scores[retained]
        current_groups = groups[retained]
        differences.append(
            float(
                current_scores[current_groups == positive].mean()
                - current_scores[current_groups == negative].mean()
            )
        )
    return bool(all(value > 0 for value in differences)), min(differences), max(differences)


def map_human_modules(
    modules: dict[str, pd.DataFrame],
    lookups,
    target_ids: set[str],
) -> dict[str, pd.DataFrame]:
    result = {}
    for module_id in MODULE_IDS:
        source = M12.add_resolution(modules[module_id], "gene", lookups)
        source = source[source["human_gene_id"].notna()].drop_duplicates("human_gene_id")
        result[module_id] = source[source["human_gene_id"].isin(target_ids)].copy()
    return result


def rank_score(expression: pd.DataFrame, mapped: pd.DataFrame, id_column: str) -> pd.Series:
    ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
    up = mapped.loc[mapped["direction"] == "up", id_column].astype(str)
    down = mapped.loc[mapped["direction"] == "down", id_column].astype(str)
    up = [gene for gene in up if gene in ranks.index]
    down = [gene for gene in down if gene in ranks.index]
    if len(up) < 10 or len(down) < 10:
        return pd.Series(np.nan, index=expression.columns)
    return ranks.loc[up].mean(axis=0) - ranks.loc[down].mean(axis=0)


def gene_matched_p(
    mapped: pd.DataFrame,
    universe: pd.DataFrame,
    id_column: str,
    effect_column: str,
    expression_column: str,
    rng: np.random.Generator,
) -> tuple[float, float]:
    universe = universe[[id_column, effect_column, expression_column]].dropna().copy()
    universe[id_column] = universe[id_column].astype(str)
    universe["expression_decile"] = pd.qcut(
        universe[expression_column].rank(method="first"), 10, labels=False
    ).astype(int)
    decile_by_id = universe.set_index(id_column)["expression_decile"]
    mapped = mapped.copy()
    mapped[id_column] = mapped[id_column].astype(str)
    mapped["expression_decile"] = mapped[id_column].map(decile_by_id)
    mapped = mapped[mapped["expression_decile"].notna()].copy()
    effect_by_id = universe.set_index(id_column)[effect_column]
    mapped["target_effect"] = mapped[id_column].map(effect_by_id)
    up = mapped[mapped["direction"] == "up"]
    down = mapped[mapped["direction"] == "down"]
    if len(up) < 10 or len(down) < 10:
        return np.nan, np.nan
    effects = {
        int(decile): block[effect_column].to_numpy(dtype=float)
        for decile, block in universe.groupby("expression_decile")
    }
    up_counts = up["expression_decile"].value_counts().to_dict()
    down_counts = down["expression_decile"].value_counts().to_dict()
    observed = float(up["target_effect"].mean() - down["target_effect"].mean())
    null = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        up_sum = sum(
            float(rng.choice(effects[int(decile)], int(count), replace=False).sum())
            for decile, count in up_counts.items()
        )
        down_sum = sum(
            float(rng.choice(effects[int(decile)], int(count), replace=False).sum())
            for decile, count in down_counts.items()
        )
        null[iteration] = up_sum / len(up) - down_sum / len(down)
    p_value = float((1 + np.sum(null >= observed - 1e-12)) / (N_NULL + 1))
    return observed, p_value


def score_matched_p(
    mapped: pd.DataFrame,
    expression: pd.DataFrame,
    mean_expression: pd.Series,
    id_column: str,
    groups: pd.Series,
    positive: str,
    negative: str,
    rng: np.random.Generator,
) -> float:
    ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
    mean_expression = mean_expression.loc[ranks.index]
    deciles = pd.qcut(mean_expression.rank(method="first"), 10, labels=False).astype(int)
    mapped = mapped.copy()
    mapped[id_column] = mapped[id_column].astype(str)
    mapped["expression_decile"] = mapped[id_column].map(deciles)
    mapped = mapped[mapped["expression_decile"].notna()]
    up = mapped[mapped["direction"] == "up"]
    down = mapped[mapped["direction"] == "down"]
    if len(up) < 10 or len(down) < 10:
        return np.nan
    observed_scores = ranks.loc[up[id_column]].mean(axis=0) - ranks.loc[down[id_column]].mean(axis=0)
    observed = float(
        observed_scores[groups == positive].mean() - observed_scores[groups == negative].mean()
    )
    ids_by_decile = {
        int(decile): block.index.to_numpy(dtype=str) for decile, block in ranks.groupby(deciles)
    }
    up_counts = up["expression_decile"].value_counts().to_dict()
    down_counts = down["expression_decile"].value_counts().to_dict()
    extreme = 0
    for _ in range(N_NULL):
        sampled_up = np.concatenate(
            [
                rng.choice(ids_by_decile[int(decile)], int(count), replace=False)
                for decile, count in up_counts.items()
            ]
        )
        sampled_down = np.concatenate(
            [
                rng.choice(ids_by_decile[int(decile)], int(count), replace=False)
                for decile, count in down_counts.items()
            ]
        )
        null_scores = ranks.loc[sampled_up].mean(axis=0) - ranks.loc[sampled_down].mean(axis=0)
        null_difference = float(
            null_scores[groups == positive].mean() - null_scores[groups == negative].mean()
        )
        extreme += int(null_difference >= observed - 1e-12)
    return float((1 + extreme) / (N_NULL + 1))


def component_counts(mapped: pd.DataFrame) -> tuple[int, int, bool]:
    up = int((mapped["direction"] == "up").sum())
    down = int((mapped["direction"] == "down").sum())
    return up, down, bool(up >= 10 and down >= 10)


def analyze_independent_human(modules, lookups):
    universe, expression, audit = M06.load_target(lookups)
    metadata, metadata_audit = M06.load_metadata(list(expression.columns))
    target_ids = set(universe["human_gene_id"].astype(str))
    mapped_modules = map_human_modules(modules, lookups, target_ids)
    universe = universe.copy()
    universe["human_gene_id"] = universe["human_gene_id"].astype(str)
    expression.index = expression.index.astype(str)
    ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
    group_series = metadata.set_index("donor").loc[expression.columns, "group"]
    rows = []
    scores_long = []
    mapping_long = []
    for index, module_id in enumerate(MODULE_IDS):
        mapped = mapped_modules[module_id].merge(
            universe[["human_gene_id", "stat", "log2FoldChange", "baseMean"]],
            on="human_gene_id",
            how="inner",
        )
        mapped["expression_decile"] = mapped["human_gene_id"].map(
            pd.qcut(universe.set_index("human_gene_id")["baseMean"].rank(method="first"), 10, labels=False)
        )
        up, down, enough = component_counts(mapped)
        if enough:
            _, gene_p = M06.expression_matched_null(
                mapped[mapped["direction"] == "up"],
                mapped[mapped["direction"] == "down"],
                universe.assign(
                    expression_decile=pd.qcut(
                        universe["baseMean"].rank(method="first"), 10, labels=False
                    ).astype(int)
                ),
                np.random.default_rng(SEED + index),
            )
            gene_concordance = float(
                mapped.loc[mapped["direction"] == "up", "stat"].mean()
                - mapped.loc[mapped["direction"] == "down", "stat"].mean()
            )
            score = ranks.loc[mapped.loc[mapped["direction"] == "up", "human_gene_id"]].mean(axis=0) - ranks.loc[
                mapped.loc[mapped["direction"] == "down", "human_gene_id"]
            ].mean(axis=0)
            score_diff = float(score[group_series == "DPN"].mean() - score[group_series == "Control"].mean())
            effect = M06.hedges_g(score, group_series)
            exact_p, assignments, _ = M06.sex_stratified_exact_permutation(score, metadata)
            loo, loo_min, loo_max = M06.leave_one_out(score, group_series)
        else:
            gene_p = gene_concordance = score_diff = effect = exact_p = loo_min = loo_max = np.nan
            assignments = 0
            loo = False
            score = pd.Series(np.nan, index=expression.columns)
        for sample, value in score.items():
            scores_long.append(
                {
                    "target_id": "independent_human_hDRG",
                    "module_id": module_id,
                    "sample_id": sample,
                    "group": group_series[sample],
                    "score": value,
                }
            )
        mapped["target_id"] = "independent_human_hDRG"
        mapping_long.append(mapped)
        rows.append(
            {
                "target_id": "independent_human_hDRG",
                "species": "human",
                "compartment": "hDRG",
                "context": "terminal_DPN_vs_control",
                "module_id": module_id,
                "mapped_up_n": up,
                "mapped_down_n": down,
                "minimum_10_per_direction": enough,
                "gene_concordance": gene_concordance,
                "gene_matched_p": gene_p,
                "score_difference": score_diff,
                "hedges_g": effect,
                "exact_p": exact_p,
                "exact_assignments": assignments,
                "loo_all_positive": loo,
                "loo_min_difference": loo_min,
                "loo_max_difference": loo_max,
                "support_rule": "human_dual_layer",
            }
        )
    tests = pd.DataFrame(rows)
    tests["gene_matched_q"] = M12.bh_adjust(tests["gene_matched_p"])
    tests["exact_q"] = M12.bh_adjust(tests["exact_p"])
    tests["component_supportive"] = (
        tests["minimum_10_per_direction"]
        & (tests["gene_concordance"] > 0)
        & (tests["gene_matched_q"] < 0.10)
        & (tests["score_difference"] > 0)
        & (tests["exact_q"] < 0.10)
        & tests["loo_all_positive"]
    )
    audit.update(metadata_audit)
    return tests, pd.DataFrame(scores_long), pd.concat(mapping_long, ignore_index=True), audit


def analyze_rat_drg(modules, human_lookups):
    counts_path = TABLES / f"GSE176017_animal_pseudobulk_raw_counts_{DATE}.tsv.gz"
    metadata_path = PHASE / "metadata" / f"GSE176017_animal_metadata_{DATE}.tsv"
    counts = pd.read_csv(counts_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t")
    samples = metadata["sample_id"].tolist()
    rat_lookup = M07.build_lookup(M07.load_gene_info(NCBI / "Rattus_norvegicus.gene_info.gz"))
    resolved = counts["gene_symbol"].map(lambda value: M07.resolve(value, rat_lookup))
    counts["target_gene_id"] = [item[0] for item in resolved]
    counts = counts[counts["target_gene_id"].notna()].copy()
    counts["row_total"] = counts[samples].sum(axis=1)
    counts = counts.sort_values(
        ["target_gene_id", "row_total", "gene_symbol"], ascending=[True, False, True]
    ).drop_duplicates("target_gene_id")
    matrix = counts.set_index("target_gene_id")[samples].astype(float)
    library_sizes = matrix.sum(axis=0)
    cpm = matrix.divide(library_sizes, axis=1) * 1_000_000
    expressed = (cpm >= 1).sum(axis=1) >= 2
    expression = np.log2(cpm.loc[expressed] + 1)
    expression.index = expression.index.astype(str)
    mean_expression = expression.mean(axis=1)
    groups = metadata.set_index("sample_id").loc[samples, "group"]
    groups.index = samples
    orthology, orthology_audit = M07.reciprocal_one_to_one_orthologues(
        NCBI / "gene_orthologs.gz", "10116"
    )
    rows = []
    scores_long = []
    mapping_long = []
    for index, module_id in enumerate(MODULE_IDS):
        source = M12.add_resolution(modules[module_id], "gene", human_lookups)
        source["target_gene_id"] = source["human_gene_id"].map(orthology)
        mapped = source[source["target_gene_id"].isin(expression.index)].drop_duplicates(
            "target_gene_id"
        )
        up, down, enough = component_counts(mapped)
        score = rank_score(expression, mapped, "target_gene_id")
        if enough:
            score_diff = float(
                score[groups == "Painful_DPN"].mean()
                - score[groups == "Diabetes_no_allodynia"].mean()
            )
            effect = hedges_g(score, groups, "Painful_DPN", "Diabetes_no_allodynia")
            exact_p, assignments = exact_permutation(
                score, groups, "Painful_DPN", "Diabetes_no_allodynia", True
            )
            loo, loo_min, loo_max = leave_one_out(
                score, groups, "Painful_DPN", "Diabetes_no_allodynia"
            )
            matched_p = score_matched_p(
                mapped,
                expression,
                mean_expression,
                "target_gene_id",
                groups,
                "Painful_DPN",
                "Diabetes_no_allodynia",
                np.random.default_rng(SEED + 200 + index),
            )
        else:
            score_diff = effect = exact_p = loo_min = loo_max = matched_p = np.nan
            assignments = 0
            loo = False
        for sample, value in score.items():
            scores_long.append(
                {
                    "target_id": "GSE176017_rat_DRG",
                    "module_id": module_id,
                    "sample_id": sample,
                    "group": groups[sample],
                    "score": value,
                }
            )
        mapped["target_id"] = "GSE176017_rat_DRG"
        mapping_long.append(mapped)
        rows.append(
            {
                "target_id": "GSE176017_rat_DRG",
                "species": "rat",
                "compartment": "DRG",
                "context": "painful_DPN_vs_diabetes_no_allodynia",
                "module_id": module_id,
                "mapped_up_n": up,
                "mapped_down_n": down,
                "minimum_10_per_direction": enough,
                "gene_concordance": np.nan,
                "gene_matched_p": np.nan,
                "score_difference": score_diff,
                "hedges_g": effect,
                "exact_p": exact_p,
                "exact_assignments": assignments,
                "loo_all_positive": loo,
                "loo_min_difference": loo_min,
                "loo_max_difference": loo_max,
                "matched_signature_p": matched_p,
                "support_rule": "rodent_progression_score",
            }
        )
    tests = pd.DataFrame(rows)
    tests["matched_signature_q"] = M12.bh_adjust(tests["matched_signature_p"])
    tests["component_supportive"] = (
        tests["minimum_10_per_direction"]
        & (tests["score_difference"] > 0)
        & (tests["hedges_g"] >= 0.5)
        & tests["loo_all_positive"]
        & (tests["matched_signature_q"] < 0.10)
    )
    audit = {
        "counts_file_sha256": M12.sha256(counts_path),
        "animal_count": len(samples),
        "expressed_gene_ids": int(len(expression)),
        "orthology": orthology_audit,
    }
    return tests, pd.DataFrame(scores_long), pd.concat(mapping_long, ignore_index=True), audit


def analyze_ocular(modules, human_lookups):
    all_tests = []
    all_scores = []
    all_mapping = []
    audits = []
    orthology_cache = {}
    target_lookup_cache = {}
    for dataset_index, (dataset, config) in enumerate(M07.DATASETS.items()):
        species = str(config["species"])
        tax_id = str(config["tax_id"])
        if species not in target_lookup_cache:
            gene_file = "Mus_musculus.gene_info.complete.gz" if species == "mouse" else "Rattus_norvegicus.gene_info.gz"
            target_lookup_cache[species] = M07.build_lookup(M07.load_gene_info(NCBI / gene_file))
            orthology_cache[species] = M07.reciprocal_one_to_one_orthologues(
                NCBI / "gene_orthologs.gz", tax_id
            )
        universe, expression, sample_meta, audit = M07.prepare_target(
            dataset, config, target_lookup_cache[species]
        )
        expression.index = expression.index.astype(str)
        groups = sample_meta.set_index("sample").loc[expression.columns, "group"]
        mean_expression = universe.set_index("target_gene_id")["mean_expression"]
        mean_expression.index = mean_expression.index.astype(str)
        orthology, orthology_audit = orthology_cache[species]
        rows = []
        for module_index, module_id in enumerate(MODULE_IDS):
            source = M12.add_resolution(modules[module_id], "gene", human_lookups)
            source["target_gene_id"] = source["human_gene_id"].map(orthology)
            mapped = source[source["target_gene_id"].isin(expression.index)].drop_duplicates(
                "target_gene_id"
            )
            up, down, enough = component_counts(mapped)
            score = rank_score(expression, mapped, "target_gene_id")
            if enough:
                score_diff = float(score[groups == "disease"].mean() - score[groups == "control"].mean())
                effect = hedges_g(score, groups, "disease", "control")
                exact_p, assignments = exact_permutation(score, groups, "disease", "control", False)
                loo, loo_min, loo_max = leave_one_out(score, groups, "disease", "control")
                gene_concordance, gene_p = gene_matched_p(
                    mapped,
                    universe,
                    "target_gene_id",
                    "log2FC",
                    "mean_expression",
                    np.random.default_rng(SEED + 500 + dataset_index * 10 + module_index),
                )
            else:
                score_diff = effect = exact_p = loo_min = loo_max = gene_concordance = gene_p = np.nan
                assignments = 0
                loo = False
            for sample, value in score.items():
                all_scores.append(
                    {
                        "target_id": dataset,
                        "module_id": module_id,
                        "sample_id": sample,
                        "group": groups[sample],
                        "score": value,
                    }
                )
            mapped["target_id"] = dataset
            all_mapping.append(mapped)
            rows.append(
                {
                    "target_id": dataset,
                    "species": species,
                    "compartment": config["compartment"],
                    "context": config["context"],
                    "module_id": module_id,
                    "mapped_up_n": up,
                    "mapped_down_n": down,
                    "minimum_10_per_direction": enough,
                    "gene_concordance": gene_concordance,
                    "gene_matched_p": gene_p,
                    "score_difference": score_diff,
                    "hedges_g": effect,
                    "exact_p": exact_p,
                    "exact_assignments": assignments,
                    "loo_all_positive": loo,
                    "loo_min_difference": loo_min,
                    "loo_max_difference": loo_max,
                    "support_rule": "ocular_gene_plus_score",
                }
            )
        current = pd.DataFrame(rows)
        current["gene_matched_q"] = M12.bh_adjust(current["gene_matched_p"])
        current["component_supportive"] = (
            current["minimum_10_per_direction"]
            & (current["score_difference"] > 0)
            & (current["hedges_g"] >= 0.5)
            & current["loo_all_positive"]
            & (current["gene_matched_q"] < 0.10)
        )
        all_tests.append(current)
        audit["orthology"] = orthology_audit
        audits.append(audit)
    return (
        pd.concat(all_tests, ignore_index=True),
        pd.DataFrame(all_scores),
        pd.concat(all_mapping, ignore_index=True),
        audits,
    )


def load_jci_reference() -> pd.DataFrame:
    gate = pd.read_csv(TABLES / f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", sep="\t")
    gate = gate[
        (gate["target_id"] == "DPN_vs_control") & gate["module_id"].isin(MODULE_IDS)
    ].copy()
    return pd.DataFrame(
        {
            "target_id": "JCI184075_human_sural_nerve",
            "species": "human",
            "compartment": "sural_nerve",
            "context": "DPN_vs_control",
            "module_id": gate["module_id"],
            "mapped_up_n": gate["mapped_up_n_gene"],
            "mapped_down_n": gate["mapped_down_n_gene"],
            "minimum_10_per_direction": gate["minimum_10_per_direction_gene"],
            "gene_concordance": gate["gene_concordance"],
            "gene_matched_p": gate["matched_null_one_sided_p"],
            "gene_matched_q": gate["matched_null_bh_q"],
            "score_difference": gate["score_mean_difference"],
            "hedges_g": gate["hedges_g"],
            "exact_p": gate["exact_one_sided_p"],
            "exact_q": gate["exact_bh_q"],
            "loo_all_positive": gate["loo_all_positive"],
            "loo_min_difference": gate["loo_min_difference"],
            "component_supportive": gate["component_pass"],
            "support_rule": "JCI_frozen_dual_layer_reference",
        }
    )


def make_heatmap(tests: pd.DataFrame, output: Path) -> None:
    target_order = [
        "independent_human_hDRG",
        "JCI184075_human_sural_nerve",
        "GSE176017_rat_DRG",
        "GSE329119",
        "GSE227165",
        "GSE180490",
    ]
    target_labels = {
        "independent_human_hDRG": "Human hDRG\nDPN vs control",
        "JCI184075_human_sural_nerve": "Human sural nerve\nDPN vs control",
        "GSE176017_rat_DRG": "Rat DRG\npainful DPN progression",
        "GSE329119": "Mouse TG\nresting diabetes",
        "GSE227165": "Rat cornea\nresting diabetes",
        "GSE180490": "Mouse cornea\nwounded diabetes",
    }
    pivot = tests.pivot(index="target_id", columns="module_id", values="hedges_g").reindex(
        index=target_order, columns=MODULE_IDS
    )
    support = tests.pivot(
        index="target_id", columns="module_id", values="component_supportive"
    ).reindex(index=target_order, columns=MODULE_IDS)
    values = pivot.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    limit = max(2.0, float(np.nanpercentile(np.abs(finite), 95))) if len(finite) else 2.0
    limit = min(limit, 3.0)
    fig, axis = plt.subplots(figsize=(9.2, 7.2))
    image = axis.imshow(
        np.clip(values, -limit, limit),
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        aspect="auto",
    )
    axis.set_xticks(range(len(MODULE_IDS)))
    axis.set_xticklabels([MODULE_DISPLAY[module] for module in MODULE_IDS], fontsize=10)
    axis.set_yticks(range(len(target_order)))
    axis.set_yticklabels([target_labels[target] for target in target_order], fontsize=10)
    for row in range(len(target_order)):
        for column in range(len(MODULE_IDS)):
            value = values[row, column]
            if not np.isfinite(value):
                label = "NA"
            else:
                label = f"g={value:.2f}"
                if bool(support.iloc[row, column]):
                    label += "\nPASS"
            color = "white" if np.isfinite(value) and abs(value) > limit * 0.55 else "black"
            axis.text(column, row, label, ha="center", va="center", fontsize=9, color=color)
    axis.set_title(
        "Source-defined hDRG components show tissue- and species-bounded transport",
        fontsize=13,
        fontweight="bold",
        pad=16,
    )
    cbar = fig.colorbar(image, ax=axis, shrink=0.82, pad=0.025)
    cbar.set_label("Hedges g (source-defined direction)")
    axis.set_xticks(np.arange(-0.5, len(MODULE_IDS), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(target_order), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=2)
    axis.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    human_lookups = M12.build_symbol_lookup(
        M12.load_gene_info(NCBI / "Homo_sapiens.gene_info.gz")
    )
    modules, _, source_audit = M12.build_source_modules()

    human_tests, human_scores, human_mapping, human_audit = analyze_independent_human(
        modules, human_lookups
    )
    rat_tests, rat_scores, rat_mapping, rat_audit = analyze_rat_drg(modules, human_lookups)
    ocular_tests, ocular_scores, ocular_mapping, ocular_audit = analyze_ocular(
        modules, human_lookups
    )
    jci_reference = load_jci_reference()

    tests = pd.concat(
        [human_tests, jci_reference, rat_tests, ocular_tests], ignore_index=True, sort=False
    )
    scores = pd.concat([human_scores, rat_scores, ocular_scores], ignore_index=True)
    mapping = pd.concat([human_mapping, rat_mapping, ocular_mapping], ignore_index=True, sort=False)

    tests.to_csv(
        TABLES / f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv",
        sep="\t",
        index=False,
    )
    scores.to_csv(
        TABLES / f"cross_target_hDRG_component_sample_scores_{DATE}.tsv",
        sep="\t",
        index=False,
    )
    mapping.to_csv(
        TABLES / f"cross_target_hDRG_component_mapping_{DATE}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    qc = {
        "analysis_date": DATE,
        "random_seed": SEED,
        "matched_null_iterations": N_NULL,
        "analysis_role": "post_JCI_explanatory_not_independent_confirmation",
        "source_audit": source_audit,
        "independent_human_hDRG_audit": human_audit,
        "rat_DRG_audit": rat_audit,
        "ocular_audits": ocular_audit,
        "support_counts_by_target": tests.groupby("target_id")["component_supportive"]
        .sum()
        .astype(int)
        .to_dict(),
    }
    with (TABLES / f"cross_target_hDRG_component_transportability_qc_{DATE}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    make_heatmap(
        tests,
        FIGURES / f"cross_target_hDRG_component_transportability_atlas_{DATE}",
    )
    print(
        tests[
            [
                "target_id",
                "module_id",
                "mapped_up_n",
                "mapped_down_n",
                "score_difference",
                "hedges_g",
                "gene_matched_q",
                "matched_signature_q",
                "exact_q",
                "loo_all_positive",
                "component_supportive",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
