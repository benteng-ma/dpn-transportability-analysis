#!/usr/bin/env python3
"""Build manuscript-ready composite Figures 1, 2, 5, and 6 from audited tables."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
TABLES = PHASE / "results" / "tables"
OUT = PHASE / "results" / "figures" / "main_composites"
OUT.mkdir(parents=True, exist_ok=True)
DATE = "2026-08-27"

BLUE = "#3C78A8"
RED = "#D75B58"
TEAL = "#2A9D8F"
ORANGE = "#E9A23B"
GREY = "#6B7280"
LIGHT_GREY = "#E5E7EB"
NAVY = "#1F3B5B"
DARK = "#17212B"
PALE_BLUE = "#E8F1F7"
PALE_RED = "#FAECEA"
PALE_TEAL = "#E8F5F2"


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, sep="\t")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.6, alpha=0.65)
    ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.06, label, transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}_{DATE}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}_{DATE}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_figure1() -> None:
    sig = read(f"hDRG_frozen_primary_stage_signatures_{DATE}.tsv")
    comp = read(f"hDRG_source_defined_transport_components_{DATE}.tsv")

    fig = plt.figure(figsize=(15.5, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # A: evidence hierarchy.
    ax_a.set_axis_off()
    panel_label(ax_a, "A")
    ax_a.set_title("Frozen transportability design", loc="left", pad=6)
    boxes = [
        (0.04, 0.60, 0.25, 0.25, "Discovery source", "Human hDRG snRNA/spatial\n2026 preprint", PALE_RED, RED),
        (0.38, 0.60, 0.25, 0.25, "Human neural targets", "Independent hDRG\nHuman sural nerve", PALE_TEAL, TEAL),
        (0.71, 0.60, 0.25, 0.25, "Explanatory atlas", "Rat DRG · Mouse TG\nResting/wounded cornea", PALE_BLUE, BLUE),
        (0.38, 0.14, 0.25, 0.25, "Boundary tests", "PBMC · Whole blood\nTear proteome", "#F3F4F6", GREY),
    ]
    for x, y, w, h, title, body, fc, ec in boxes:
        ax_a.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012", fc=fc, ec=ec, lw=1.8))
        ax_a.text(x + w / 2, y + h * 0.68, title, ha="center", va="center", fontweight="bold", color=DARK)
        ax_a.text(x + w / 2, y + h * 0.32, body, ha="center", va="center", fontsize=9, color=DARK)
    arrow = dict(arrowstyle="-|>", lw=1.8, color=NAVY, mutation_scale=13)
    ax_a.annotate("", xy=(0.37, 0.72), xytext=(0.29, 0.72), arrowprops=arrow)
    ax_a.annotate("", xy=(0.70, 0.72), xytext=(0.63, 0.72), arrowprops=arrow)
    ax_a.annotate("", xy=(0.50, 0.40), xytext=(0.50, 0.59), arrowprops=arrow)
    ax_a.text(0.50, 0.04, "Membership, direction and gates frozen before target scoring", ha="center", fontsize=9.5, fontweight="bold", color=NAVY)
    ax_a.text(0.04, 0.92, "Evidence order", fontsize=9, color=GREY)

    # B: five frozen signatures.
    panel_label(ax_b, "B")
    labels = {
        "early_allcell_diabetes_vs_control": "Early all-cell\nDM vs control",
        "late_allcell_DPN_vs_diabetes": "Late all-cell\nDPN vs DM",
        "late_neuron_DPN_vs_diabetes": "Late neuronal\nDPN vs DM",
        "severity_neuron_modhigh_vs_low_nageotte": "Neuronal severity\nmod/high vs low",
        "xenium_DPN_vs_control": "Spatial Xenium\nDPN vs control",
    }
    order = list(labels)
    counts = sig.groupby(["contrast_id", "direction"]).size().unstack(fill_value=0).reindex(order)
    y = np.arange(len(order))
    ax_b.barh(y, counts["up"], color=RED, label="Source-up")
    ax_b.barh(y, counts["down"], left=counts["up"], color=BLUE, label="Source-down")
    for i, key in enumerate(order):
        total = int(counts.loc[key].sum())
        ax_b.text(total + max(counts.sum(axis=1)) * 0.018, i, f"n={total:,}", va="center", fontsize=9)
    ax_b.set_yticks(y, [labels[key] for key in order])
    ax_b.invert_yaxis()
    ax_b.set_xlabel("Frozen source genes")
    ax_b.set_title("Five directional hDRG programs", loc="left")
    ax_b.legend(frameon=False, ncol=2, loc="lower right")
    clean_axis(ax_b)

    def component_panel(ax: plt.Axes, title: str, module_order: list[str], display: dict[str, str], opposed: str) -> None:
        panel_label(ax, "C" if "Late" in title else "D")
        rows = []
        for module in module_order:
            subset = comp[comp["module_id"] == module]
            rows.append((display[module], int((subset["direction"] == "up").sum()), int((subset["direction"] == "down").sum())))
        opposed_n = int((comp["module_id"] == opposed).sum())
        names = [r[0] for r in rows] + ["Opposed overlap\n(audit only)"]
        ups = [r[1] for r in rows] + [opposed_n]
        downs = [r[2] for r in rows] + [0]
        yy = np.arange(len(names))
        ax.barh(yy, ups, color=[RED] * len(rows) + ["white"], edgecolor=[RED] * len(rows) + [GREY], hatch=[None] * len(rows) + ["///"], label="Up")
        ax.barh(yy, downs, left=ups, color=BLUE, label="Down")
        for i, (u, d) in enumerate(zip(ups, downs)):
            ax.text(u + d + max(np.array(ups) + np.array(downs)) * 0.02, i, f"{u+d:,}", va="center", fontsize=9)
        ax.set_yticks(yy, names)
        ax.invert_yaxis()
        ax.set_xlabel("Source-defined genes")
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, ncol=2, loc="lower right")
        clean_axis(ax)

    component_panel(
        ax_c,
        "Late-stage source-only decomposition",
        ["late_shared_concordant_neuronal_core", "late_neuron_residual", "late_allcell_residual"],
        {
            "late_shared_concordant_neuronal_core": "Shared neuronal core",
            "late_neuron_residual": "Neuron residual",
            "late_allcell_residual": "All-cell residual",
        },
        "late_directionally_opposed_overlap",
    )
    component_panel(
        ax_d,
        "Severity source-only decomposition",
        ["severity_neuron_shared_concordant_core", "severity_neuron_residual"],
        {
            "severity_neuron_shared_concordant_core": "Severity-shared core",
            "severity_neuron_residual": "Severity residual",
        },
        "severity_directionally_opposed_overlap",
    )
    fig.suptitle("Human hDRG disease programs were frozen and decomposed before transport testing", fontsize=17, fontweight="bold")
    save(fig, "Figure1_source_design_and_decomposition")


def build_figure2() -> None:
    tests = read(f"independent_human_DPN_bulk_signature_tests_{DATE}.tsv")
    scores = read(f"independent_human_DPN_bulk_donor_scores_{DATE}.tsv")
    cross = read(f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv")
    cells = read(f"hDRG_major_celltype_pairwise_proportion_tests_{DATE}.tsv")
    severity_id = "severity_neuron_modhigh_vs_low_nageotte"
    row = tests.loc[tests["contrast_id"] == severity_id].iloc[0]
    sev_scores = scores.loc[scores["contrast_id"] == severity_id].copy()

    fig = plt.figure(figsize=(15.5, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_label(ax_a, "A")
    qvals = np.array([row["gene_matched_bh_q"], row["donor_exact_bh_q"]], dtype=float)
    vals = -np.log10(qvals)
    bars = ax_a.bar([0, 1], vals, color=[TEAL, ORANGE], width=0.58)
    ax_a.axhline(-np.log10(0.10), color=GREY, ls="--", lw=1.2, label="Frozen Q=0.10 gate")
    ax_a.set_xticks([0, 1], ["Gene-level\nmatched null", "Donor-level\nsex-stratified exact"])
    ax_a.set_ylabel("-log10(Q)")
    ax_a.set_title("Original neuronal-severity program passes both layers", loc="left")
    for bar, q in zip(bars, qvals):
        ax_a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.08, f"Q={q:.4g}", ha="center", fontsize=9, fontweight="bold")
    ax_a.text(0.98, 0.96, f"Mapped: {int(row['mapped_up_n'])} up / {int(row['mapped_down_n'])} down\nHedges g={row['donor_hedges_g']:.3f}\nAll leave-one-out positive", transform=ax_a.transAxes, ha="right", va="top", bbox=dict(boxstyle="round", fc=PALE_TEAL, ec=TEAL), fontsize=9)
    clean_axis(ax_a)
    ax_a.set_ylim(0, max(vals) + 0.25)
    ax_a.legend(frameon=False, loc="lower right")

    panel_label(ax_b, "B")
    rng = np.random.default_rng(20260827)
    for x, group, color in [(0, "Control", BLUE), (1, "DPN", RED)]:
        subset = sev_scores[sev_scores["group"] == group]
        jitter = rng.uniform(-0.08, 0.08, size=len(subset))
        for j, (_, r) in enumerate(subset.iterrows()):
            marker = "^" if r["Sex"] == "M" else "o"
            ax_b.scatter(x + jitter[j], r["score"], s=58, marker=marker, color=color, edgecolor="white", linewidth=0.6, zorder=3)
        mean = subset["score"].mean()
        sem = subset["score"].sem()
        ax_b.errorbar(x, mean, yerr=sem, fmt="_", markersize=20, capsize=5, color=DARK, lw=1.8, zorder=4)
    ax_b.set_xticks([0, 1], ["Control\n(n=7)", "Painful DPN\n(n=5)"])
    ax_b.set_ylabel("Within-donor rank score")
    ax_b.set_title("Independent human hDRG donor scores", loc="left")
    ax_b.text(0.03, 0.95, f"Δ={row['donor_score_difference']:.5f}\nexact Q={row['donor_exact_bh_q']:.4f}", transform=ax_b.transAxes, va="top", fontsize=9)
    ax_b.text(0.97, 0.05, "○ female   △ male", transform=ax_b.transAxes, ha="right", fontsize=8.5, color=GREY)
    clean_axis(ax_b)

    panel_label(ax_c, "C")
    h = cross[cross["target_id"] == "independent_human_hDRG"].copy()
    order = ["late_shared_concordant_neuronal_core", "late_neuron_residual", "late_allcell_residual"]
    labels = ["Shared neuronal core", "Neuron residual", "All-cell residual"]
    h = h.set_index("module_id").loc[order]
    g = h["hedges_g"].astype(float).to_numpy()
    supportive = h["component_supportive"].astype(str).str.lower().eq("true").to_numpy()
    colors = np.where(supportive, TEAL, GREY)
    yy = np.arange(3)
    ax_c.barh(yy, g, color=colors)
    ax_c.axvline(0, color=DARK, lw=1)
    ax_c.set_yticks(yy, labels)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("Hedges g (DPN - control)")
    ax_c.set_title("Post-primary component comparison", loc="left")
    for i, (val, ok, qg, qe) in enumerate(zip(g, supportive, h["gene_matched_q"], h["exact_q"])):
        text_x = val + 0.04 if val >= 0 else 0.03
        ax_c.text(text_x, i, f"g={val:.2f}  {'PASS' if ok else 'FAIL'}\ngene Q={float(qg):.3g}; exact Q={float(qe):.3g}", va="center", ha="left", fontsize=8.4)
    ax_c.set_xlim(-0.48, 1.32)
    clean_axis(ax_c)

    panel_label(ax_d, "D")
    ax_d.set_axis_off()
    ax_d.set_title("Source hDRG cell-proportion context: diabetes to DPN", loc="left")
    selected = cells[(cells["group1"] == "Diabetic") & (cells["group2"] == "DPN") & (cells["pairwise_q_below_0_05"].astype(str).str.lower() == "true")].copy()
    selected["display"] = selected["cell_type"].replace({"SMC": "Smooth-muscle cells"})
    y_positions = [0.72, 0.47, 0.22]
    for (_, r), y in zip(selected.iterrows(), y_positions):
        up = "increased" in str(r["source_article_narrative_for_DPN"])
        color = RED if up else BLUE
        symbol = "↑" if up else "↓"
        ax_d.add_patch(patches.FancyBboxPatch((0.08, y - 0.09), 0.80, 0.16, boxstyle="round,pad=0.012", fc=PALE_RED if up else PALE_BLUE, ec=color, lw=1.5))
        ax_d.text(0.16, y, symbol, fontsize=24, color=color, va="center", ha="center", fontweight="bold")
        ax_d.text(0.26, y + 0.025, r["display"], fontsize=11, fontweight="bold", va="center")
        ax_d.text(0.26, y - 0.035, f"source pairwise Q={float(r['p.adj']):.4g}", fontsize=9, color=GREY, va="center")
    ax_d.text(0.08, 0.03, "Context only: these source proportions do not identify the cell origin of a bulk target score.", fontsize=9, color=GREY)

    fig.suptitle("A terminal human hDRG degeneration program reproduces in an independent donor cohort", fontsize=17, fontweight="bold")
    save(fig, "Figure2_independent_human_hDRG_validation")


def build_figure5() -> None:
    cross = read(f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv")
    order_targets = ["independent_human_hDRG", "JCI184075_human_sural_nerve", "GSE176017_rat_DRG", "GSE329119", "GSE227165", "GSE180490"]
    order_modules = ["late_shared_concordant_neuronal_core", "late_neuron_residual", "late_allcell_residual"]
    target_labels = ["Human hDRG\nDPN vs control", "Human sural nerve\nDPN vs control", "Rat DRG\npainful progression", "Mouse TG\nresting diabetes", "Rat cornea\nresting diabetes", "Mouse cornea\nwounded diabetes"]
    module_labels = ["Shared neuronal core", "Neuron residual", "All-cell residual"]
    pivot = cross.pivot(index="target_id", columns="module_id", values="hedges_g").loc[order_targets, order_modules].astype(float)
    status = cross.set_index(["target_id", "module_id"])["component_supportive"]

    fig = plt.figure(figsize=(15.5, 9.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.1, 1.1, 0.9], height_ratios=[1.1, 1.0])
    ax_a = fig.add_subplot(gs[:, :2])
    ax_b = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[1, 2])

    panel_label(ax_a, "A")
    norm = TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
    im = ax_a.imshow(pivot.to_numpy(), cmap="RdBu_r", norm=norm, aspect="auto")
    ax_a.set_xticks(np.arange(3), module_labels)
    ax_a.set_yticks(np.arange(6), target_labels)
    ax_a.set_title("Component effect atlas", loc="left")
    for i, target in enumerate(order_targets):
        for j, module in enumerate(order_modules):
            value = float(pivot.loc[target, module])
            ok = str(status.loc[(target, module)]).lower() == "true"
            tier = "PASS" if target.startswith(("independent", "JCI")) else "provisional"
            text = f"g={value:.2f}" + (f"\n{tier}" if ok else "")
            ax_a.text(j, i, text, ha="center", va="center", color="white" if abs(value) > 1.0 else DARK, fontsize=9, fontweight="bold" if ok else "normal")
    for edge in np.arange(-0.5, 6, 1):
        ax_a.axhline(edge, color="white", lw=2)
    for edge in np.arange(-0.5, 3, 1):
        ax_a.axvline(edge, color="white", lw=2)
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.035, pad=0.03)
    cbar.set_label("Hedges g in source-defined direction")

    panel_label(ax_b, "B")
    ax_b.set_axis_off()
    ax_b.set_title("Evidence key", loc="left")
    key_rows = [
        (0.78, TEAL, "Human PASS", "Frozen gene + biological-unit + LOO gate"),
        (0.52, ORANGE, "Provisional", "Post-primary animal/ocular support"),
        (0.26, GREY, "No support", "Failed, reversed, or insufficient gate"),
    ]
    for y, color, head, detail in key_rows:
        ax_b.add_patch(patches.FancyBboxPatch((0.04, y - 0.08), 0.90, 0.16, boxstyle="round,pad=0.01", fc="white", ec=color, lw=1.8))
        ax_b.add_patch(patches.Rectangle((0.07, y - 0.035), 0.06, 0.07, fc=color, ec="none"))
        ax_b.text(0.17, y + 0.025, head, fontweight="bold", va="center")
        ax_b.text(0.17, y - 0.030, detail, fontsize=8.5, color=GREY, va="center")
    ax_b.text(0.04, 0.04, "Human targets carry the confirmatory weight.\nAnimal and ocular rows explain boundaries.", fontsize=9, color=NAVY, fontweight="bold")

    panel_label(ax_c, "C")
    cornea = pivot.loc[["GSE227165", "GSE180490"]]
    xx = np.arange(3)
    width = 0.35
    ax_c.bar(xx - width / 2, cornea.iloc[0], width, color=BLUE, label="Resting diabetic cornea")
    ax_c.bar(xx + width / 2, cornea.iloc[1], width, color=RED, label="Wounded diabetic cornea")
    ax_c.axhline(0, color=DARK, lw=1)
    ax_c.set_xticks(xx, ["Shared\ncore", "Neuron\nresidual", "All-cell\nresidual"])
    ax_c.set_ylabel("Hedges g")
    ax_c.set_title("Corneal injury-context boundary", loc="left")
    ax_c.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax_c.text(1.0 + width / 2, cornea.iloc[1, 1] - 0.08, "provisional", ha="center", color="white", fontsize=8.0, fontweight="bold", rotation=90, va="top")
    ax_c.text(2.0 + width / 2, cornea.iloc[1, 2] - 0.08, "provisional", ha="center", color="white", fontsize=8.0, fontweight="bold", rotation=90, va="top")
    clean_axis(ax_c)

    fig.suptitle("Source-defined hDRG components show anatomy-, species-, and injury-bounded transport", fontsize=17, fontweight="bold")
    save(fig, "Figure5_cross_target_transportability_atlas")


def build_figure6() -> None:
    pbmc = read(f"human_PBMC_stage_projection_meta_analysis_{DATE}.tsv")
    clinical = read(f"GSE302658_clinical_signature_tests_{DATE}.tsv")
    tear = read(f"PXD062366_hDRG_signature_coverage_{DATE}.tsv")

    fig = plt.figure(figsize=(15.5, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_label(ax_a, "A")
    labels = ["Early: diabetes vs healthy", "Late: DPN vs diabetes"]
    est = pbmc["fixed_effect_hedges_g"].astype(float).to_numpy()
    se = pbmc["standard_error"].astype(float).to_numpy()
    yy = np.arange(2)
    ax_a.errorbar(est, yy, xerr=1.96 * se, fmt="o", color=GREY, ecolor=GREY, capsize=4, markersize=8)
    ax_a.axvline(0, color=DARK, lw=1)
    ax_a.axvline(0.5, color=TEAL, lw=1, ls="--", label="Prespecified positive-effect target")
    ax_a.set_yticks(yy, labels)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Fixed-effect Hedges g (95% CI)")
    ax_a.set_title("Two-cohort PBMC stage projection", loc="left")
    for i, (g, i2) in enumerate(zip(est, pbmc["I2_percent"].astype(float))):
        y_text = i + 0.16 if i == 0 else i - 0.16
        ax_a.text(0.98, y_text, f"g={g:.2f}; I²={i2:.1f}%  FAIL", transform=ax_a.get_yaxis_transform(), ha="right", va="center", fontsize=8.5, color=GREY, fontweight="bold")
    clean_axis(ax_a)
    ax_a.text(0.5, 0.06, "Prespecified positive-effect target", transform=ax_a.get_xaxis_transform(), ha="center", fontsize=8.2, color=TEAL)

    panel_label(ax_b, "B")
    cl = clinical[clinical["test_family"].isin(["primary", "longitudinal_primary"])].copy()
    labels_b = ["Baseline NPSI", "Change in NPSI"]
    rho = cl["rho"].astype(float).to_numpy()
    lo = cl["bootstrap_ci_low"].astype(float).to_numpy()
    hi = cl["bootstrap_ci_high"].astype(float).to_numpy()
    xerr = np.vstack([rho - lo, hi - rho])
    yy = np.arange(2)
    ax_b.errorbar(rho, yy, xerr=xerr, fmt="o", color=ORANGE, ecolor=ORANGE, capsize=4, markersize=8)
    ax_b.axvline(0, color=DARK, lw=1)
    ax_b.set_yticks(yy, labels_b)
    ax_b.invert_yaxis()
    ax_b.set_xlabel("Spearman rho (bootstrap 95% CI)")
    ax_b.set_title("Clinical whole-blood transfer", loc="left")
    for i, r in enumerate(cl.itertuples(index=False)):
        matched = "n/a" if pd.isna(r.matched_signature_p) else f"{float(r.matched_signature_p):.3f}"
        y_text = i + 0.16 if i == 0 else i - 0.16
        ax_b.text(0.98, y_text, f"permutation P={float(r.permutation_positive_p):.3f}; adjusted P={float(r.hc3_exposure_p):.3f}; matched P={matched} → FAIL", transform=ax_b.get_yaxis_transform(), ha="right", va="center", fontsize=8.1, color=GREY)
    clean_axis(ax_b)

    panel_label(ax_c, "C")
    program_order = [
        "early_allcell_diabetes_vs_control",
        "late_allcell_DPN_vs_diabetes",
        "late_neuron_DPN_vs_diabetes",
        "severity_neuron_modhigh_vs_low_nageotte",
    ]
    labels_c = ["Early all-cell", "Late all-cell", "Late neuronal", "Neuronal severity"]
    t = tear[tear["contrast_id"].isin(program_order)].pivot(index="contrast_id", columns="direction", values="tear_proteome_gene_n_detection_ge_50pct").reindex(program_order)
    xx = np.arange(4)
    width = 0.36
    ax_c.bar(xx - width / 2, t["up"], width, color=RED, label="Source-up detected ≥50%")
    ax_c.bar(xx + width / 2, t["down"], width, color=BLUE, label="Source-down detected ≥50%")
    ax_c.set_xticks(xx, labels_c, rotation=15, ha="right")
    ax_c.set_ylabel("Mapped tear proteins")
    ax_c.set_title("Tear-proteome coverage stop", loc="left")
    ax_c.legend(frameon=False, fontsize=8.5)
    ax_c.text(0.98, 0.95, "Direct testing stopped:\nsparse/asymmetric directions", transform=ax_c.transAxes, ha="right", va="top", bbox=dict(boxstyle="round", fc=PALE_RED, ec=RED), fontsize=9, fontweight="bold")
    clean_axis(ax_c)

    panel_label(ax_d, "D")
    ax_d.set_axis_off()
    ax_d.set_title("Interpretive boundary", loc="left")
    ax_d.add_patch(patches.FancyBboxPatch((0.04, 0.58), 0.38, 0.27, boxstyle="round,pad=0.015", fc=PALE_TEAL, ec=TEAL, lw=2))
    ax_d.text(0.23, 0.76, "Neural-tissue\ntransportability", ha="center", va="center", fontsize=13, fontweight="bold", color=NAVY)
    ax_d.text(0.23, 0.63, "Selective human hDRG ↔ sural\ncomponent retention", ha="center", va="center", fontsize=9)
    ax_d.add_patch(patches.FancyBboxPatch((0.58, 0.58), 0.38, 0.27, boxstyle="round,pad=0.015", fc="#F3F4F6", ec=GREY, lw=2))
    ax_d.text(0.77, 0.76, "Accessible biomarker\nportability", ha="center", va="center", fontsize=13, fontweight="bold", color=DARK)
    ax_d.text(0.77, 0.63, "PBMC FAIL · Whole blood FAIL\nTears NOT TESTABLE", ha="center", va="center", fontsize=9)
    ax_d.text(0.50, 0.73, "≠", ha="center", va="center", fontsize=24, fontweight="bold", color=GREY)
    ax_d.text(0.50, 0.38, "A failed frozen projection rejects a direct proxy;\nit does not imply absence of all diabetes-related biology.", ha="center", va="center", fontsize=10.5, color=NAVY, fontweight="bold", bbox=dict(boxstyle="round", fc="white", ec=LIGHT_GREY))
    ax_d.text(0.50, 0.16, "Transportability and biomarker portability are separate hypotheses.", ha="center", va="center", fontsize=11.5, fontweight="bold", color=DARK)

    fig.suptitle("Accessible compartments do not provide a stable direct readout of frozen hDRG programs", fontsize=17, fontweight="bold")
    save(fig, "Figure6_accessible_compartment_boundary")


def main() -> None:
    build_figure1()
    build_figure2()
    build_figure5()
    build_figure6()
    source_figures = PHASE / "results" / "figures"
    for source_stem, target_stem in [
        ("JCI184075_hDRG_component_transportability", "Figure3_human_sural_nerve_transport_and_severity"),
        ("hDRG_component_functional_annotation", "Figure4_component_functional_annotation"),
    ]:
        for suffix in ("png", "pdf"):
            shutil.copy2(source_figures / f"{source_stem}_{DATE}.{suffix}", OUT / f"{target_stem}_{DATE}.{suffix}")
    for path in sorted(OUT.glob(f"*_{DATE}.*")):
        print(f"{path.name}\t{path.stat().st_size}")


if __name__ == "__main__":
    main()
