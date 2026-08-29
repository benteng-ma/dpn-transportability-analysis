#!/usr/bin/env python3
"""Annotate frozen hDRG transport components with competitive pathway enrichment."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


PHASE = Path(__file__).resolve().parents[2]
SCRIPTS = PHASE / "analysis" / "scripts"
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
REFERENCE = PHASE / "data" / "raw" / "reference"
DATE = "2026-08-27"
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
LIBRARIES = {
    "GO_BP_2025": REFERENCE / "GO_Biological_Process_2025.gmt",
    "Reactome_2024": REFERENCE / "Reactome_Pathways_2024.gmt",
    "Hallmark_2020": REFERENCE / "MSigDB_Hallmark_2020.gmt",
}


def load_script():
    path = SCRIPTS / "12_validate_hdrg_components_in_human_sural_nerve.py"
    spec = importlib.util.spec_from_file_location("phase06_sural", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


M12 = load_script()


def load_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            genes = {gene.strip().upper() for gene in parts[2:] if gene.strip()}
            if term and genes:
                result[term] = genes
    return result


def odds_ratio(k: int, n: int, K: int, M: int) -> float:
    a = k + 0.5
    b = n - k + 0.5
    c = K - k + 0.5
    d = M - K - n + k + 0.5
    return float((a * d) / (b * c))


def reduce_redundancy(frame: pd.DataFrame, maximum: int = 5) -> pd.DataFrame:
    retained = []
    retained_sets = []
    ordered = frame.sort_values(
        ["bh_q", "p_value", "overlap_n", "term"],
        ascending=[True, True, False, True],
    )
    for _, row in ordered.iterrows():
        genes = set(str(row["overlap_genes"]).split(";"))
        if any(
            len(genes & previous) / len(genes | previous) >= 0.50
            for previous in retained_sets
        ):
            continue
        retained.append(row)
        retained_sets.append(genes)
        if len(retained) >= maximum:
            break
    return pd.DataFrame(retained)


def wrap_term(term: str, maximum: int = 42) -> str:
    words = term.replace("_", " ").split()
    lines = []
    current = []
    for word in words:
        proposed = " ".join(current + [word])
        if current and len(proposed) > maximum:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)


def make_figure(reduced: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(14.5, 14.5))
    directions = ["up", "down"]
    for row_index, module_id in enumerate(MODULE_IDS):
        for column_index, direction in enumerate(directions):
            axis = axes[row_index, column_index]
            block = reduced[
                (reduced["module_id"] == module_id) & (reduced["direction"] == direction)
            ].copy()
            block = block.sort_values(["bh_q", "odds_ratio"]).head(7).iloc[::-1]
            if block.empty:
                axis.text(0.5, 0.5, "No Q < 0.05 term", ha="center", va="center")
                axis.set_axis_off()
                continue
            x = -np.log10(block["bh_q"].clip(lower=1e-300))
            sizes = 30 + 210 * block["overlap_n"] / block["overlap_n"].max()
            scatter = axis.scatter(
                x,
                np.arange(len(block)),
                s=sizes,
                c=np.log2(block["odds_ratio"].clip(lower=1e-6)),
                cmap="viridis",
                edgecolor="black",
                linewidth=0.4,
            )
            axis.set_yticks(np.arange(len(block)))
            axis.set_yticklabels(
                [f"{wrap_term(term)}\n[{library}]" for term, library in zip(block["term"], block["library"])],
                fontsize=8,
            )
            axis.axvline(-math.log10(0.05), color="#444444", linestyle="--", linewidth=0.8)
            axis.set_xlabel("-log10(BH Q)")
            axis.grid(axis="x", alpha=0.2)
            axis.set_title(
                f"{MODULE_DISPLAY[module_id]} — source {direction}",
                loc="left",
                fontweight="bold",
            )
            colorbar = fig.colorbar(scatter, ax=axis, fraction=0.045, pad=0.02)
            colorbar.set_label("log2 odds ratio", fontsize=8)
    fig.suptitle(
        "Competitive pathway composition of frozen hDRG transport components",
        fontsize=15,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.008,
        "Background: union of late all-cell and late-neuron source DE genes; dot size = overlap count.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    modules, _, source_audit = M12.build_source_modules()
    background = {
        gene.upper()
        for gene in pd.concat(
            [modules["original_late_allcell"], modules["original_late_neuron"]],
            ignore_index=True,
        )["gene"].astype(str)
    }
    M = len(background)
    rows = []
    library_audit = {}
    for library, path in LIBRARIES.items():
        terms = load_gmt(path)
        library_audit[library] = {
            "path": str(path),
            "sha256": M12.sha256(path),
            "deposited_terms": len(terms),
        }
        for module_id in MODULE_IDS:
            module = modules[module_id].copy()
            module["gene_upper"] = module["gene"].astype(str).str.upper()
            for direction in ["up", "down"]:
                component = set(module.loc[module["direction"] == direction, "gene_upper"]) & background
                n = len(component)
                for term, genes in terms.items():
                    term_background = genes & background
                    K = len(term_background)
                    if K < 5 or K > 500:
                        continue
                    overlap = component & term_background
                    k = len(overlap)
                    if k < 3:
                        continue
                    p_value = float(hypergeom.sf(k - 1, M, K, n))
                    rows.append(
                        {
                            "library": library,
                            "module_id": module_id,
                            "module_display": MODULE_DISPLAY[module_id],
                            "direction": direction,
                            "term": term,
                            "background_n": M,
                            "component_direction_n": n,
                            "term_background_n": K,
                            "overlap_n": k,
                            "gene_ratio": k / n,
                            "background_ratio": K / M,
                            "odds_ratio": odds_ratio(k, n, K, M),
                            "p_value": p_value,
                            "overlap_genes": ";".join(sorted(overlap)),
                        }
                    )
    results = pd.DataFrame(rows)
    results["bh_q"] = results.groupby(
        ["library", "module_id", "direction"], group_keys=False
    )["p_value"].apply(M12.bh_adjust)
    results["formal_enrichment"] = (results["bh_q"] < 0.05) & (results["odds_ratio"] > 1)
    formal = results[results["formal_enrichment"]].copy()
    top15 = formal.sort_values(
        ["library", "module_id", "direction", "bh_q", "p_value", "overlap_n"],
        ascending=[True, True, True, True, True, False],
    ).groupby(["library", "module_id", "direction"], group_keys=False).head(15)
    reduced_blocks = []
    for _, block in formal.groupby(["library", "module_id", "direction"], sort=False):
        reduced_blocks.append(reduce_redundancy(block, maximum=5))
    reduced = pd.concat(reduced_blocks, ignore_index=True) if reduced_blocks else formal.head(0)

    results.to_csv(
        TABLES / f"hDRG_component_functional_annotation_all_terms_{DATE}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    top15.to_csv(
        TABLES / f"hDRG_component_functional_annotation_top15_{DATE}.tsv",
        sep="\t",
        index=False,
    )
    reduced.to_csv(
        TABLES / f"hDRG_component_functional_annotation_redundancy_reduced_{DATE}.tsv",
        sep="\t",
        index=False,
    )
    qc = {
        "analysis_date": DATE,
        "source_audit": source_audit,
        "background_unique_genes": M,
        "library_audit": library_audit,
        "tested_term_rows": int(len(results)),
        "formal_enrichment_rows": int(len(formal)),
        "formal_counts": formal.groupby(["module_id", "direction"])
        .size()
        .to_dict(),
    }
    qc["formal_counts"] = {"|".join(key): int(value) for key, value in qc["formal_counts"].items()}
    with (TABLES / f"hDRG_component_functional_annotation_qc_{DATE}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
    make_figure(
        reduced,
        FIGURES / f"hDRG_component_functional_annotation_{DATE}",
    )
    print(qc)
    print(
        reduced[
            ["library", "module_display", "direction", "term", "overlap_n", "odds_ratio", "bh_q"]
        ].sort_values(["module_display", "direction", "bh_q"]).to_string(index=False)
    )


if __name__ == "__main__":
    main()
