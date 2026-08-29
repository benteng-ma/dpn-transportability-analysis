#!/usr/bin/env python3
"""Build the phase 0.6 manuscript figure map and numerical-claim registry."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
MANUSCRIPT = PHASE / "manuscript"
DATE = "2026-08-27"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, sep="\t")


def one(frame: pd.DataFrame, **selector: object) -> pd.Series:
    selected = frame.copy()
    for column, value in selector.items():
        selected = selected[selected[column] == value]
    assert len(selected) == 1, (selector, len(selected))
    return selected.iloc[0]


def truth(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() == "true"


def compact(values: dict[str, object]) -> str:
    cleaned = {}
    for key, value in values.items():
        if pd.isna(value):
            cleaned[key] = None
        elif isinstance(value, float):
            cleaned[key] = float(value)
        elif isinstance(value, (int, bool)):
            cleaned[key] = value
        else:
            cleaned[key] = str(value)
    return json.dumps(cleaned, ensure_ascii=False, separators=(",", ":"))


def build_figure_map() -> pd.DataFrame:
    rows = [
        {
            "figure_id": "Figure 1",
            "result_section": "Source programs and source-only decomposition",
            "panels": "A, study/evidence design; B, five frozen signatures; C, late-component decomposition; D, severity decomposition",
            "evidence_ids": "E01",
            "datasets_or_resources": "2026 human hDRG snRNA-seq/spatial preprint",
            "primary_message": "Human hDRG stage programs were frozen before target testing and decomposed without using validation cohorts.",
            "source_data": "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv; hDRG_source_defined_transport_components_2026-08-27.tsv; phase0_6_evidence_atlas_2026-08-27.tsv",
            "display_files": "new composite required",
            "evidence_label": "DISCOVERY SOURCE; NON-PEER-REVIEWED PREPRINT",
            "mandatory_caveat": "No target result was used to define membership or direction.",
            "completion_status": "SOURCE DATA READY; COMPOSITE PENDING",
        },
        {
            "figure_id": "Figure 2",
            "result_section": "Independent proximal human hDRG validation",
            "panels": "A, original severity dual-layer gate; B, donor scores/LOO; C, explanatory three-component comparison; D, source cell-proportion context",
            "evidence_ids": "E02-E05; E22",
            "datasets_or_resources": "Independent 2022 human hDRG bulk; source hDRG Sup_Data1",
            "primary_message": "Terminal degeneration severity reproduced in an independent human ganglion cohort, while late components showed selective proximal transfer.",
            "source_data": "independent_human_DPN_bulk_signature_tests_2026-08-27.tsv; independent_human_DPN_bulk_donor_scores_2026-08-27.tsv; cross_target_hDRG_component_transportability_tests_2026-08-27.tsv; hDRG_major_celltype_pairwise_proportion_tests_2026-08-27.tsv",
            "display_files": "human_DPN_bulk_signature_validation_2026-08-27.png; new explanatory/context panels required",
            "evidence_label": "INDEPENDENT HUMAN VALIDATION plus POST-PRIMARY EXPLANATORY ANALYSIS",
            "mandatory_caveat": "The validation cohort lacks a diabetes-without-neuropathy group and has age imbalance.",
            "completion_status": "CORE PANEL READY; COMPOSITE PENDING",
        },
        {
            "figure_id": "Figure 3",
            "result_section": "Selective distal human sural-nerve transfer and axonal-loss severity",
            "panels": "A, frozen gates; B, DPN-control components; C, donor scores and LOO; D, TPM sensitivity; E, severe-moderate axonal-loss components",
            "evidence_ids": "E06-E10",
            "datasets_or_resources": "JCI 2025 human sural nerve",
            "primary_message": "The late neuron residual transferred selectively to distal nerve, and a severity-shared neuronal core tracked severe axonal loss within DPN.",
            "source_data": "JCI184075_hDRG_component_gate_table_2026-08-27.tsv; JCI184075_hDRG_component_gene_tests_2026-08-27.tsv; JCI184075_hDRG_component_sample_tests_2026-08-27.tsv; JCI184075_hDRG_component_sample_scores_2026-08-27.tsv",
            "display_files": "JCI184075_hDRG_component_transportability_2026-08-27.png/pdf",
            "evidence_label": "FROZEN INDEPENDENT PEER-REVIEWED HUMAN TARGET",
            "mandatory_caveat": "Bulk sural nerve is not a purified neuronal compartment; association does not prove soma-to-axon propagation.",
            "completion_status": "READY",
        },
        {
            "figure_id": "Figure 4",
            "result_section": "Functional identities of source-defined components",
            "panels": "A, shared core; B, neuron residual; C, all-cell residual; D, interpretation model",
            "evidence_ids": "E01; E22",
            "datasets_or_resources": "Source hDRG component genes; GO BP 2025; Reactome 2024; Hallmark 2020",
            "primary_message": "The components separate sensory-neuronal signalling, neuron-intrinsic maintenance/adaptation, and broad tissue remodelling.",
            "source_data": "hDRG_component_functional_annotation_all_terms_2026-08-27.tsv.gz; hDRG_component_functional_annotation_top15_2026-08-27.tsv; hDRG_component_functional_annotation_redundancy_reduced_2026-08-27.tsv",
            "display_files": "hDRG_component_functional_annotation_2026-08-27.png/pdf",
            "evidence_label": "SOURCE-COMPONENT COMPETITIVE ANNOTATION",
            "mandatory_caveat": "Enrichment describes source composition and is not target-specific pathway activation or causal mechanism.",
            "completion_status": "READY",
        },
        {
            "figure_id": "Figure 5",
            "result_section": "Anatomical, species and injury-context transportability atlas",
            "panels": "A, Hedges-g atlas; B, formal/provisional evidence key; C, resting-versus-wounded corneal boundary",
            "evidence_ids": "E03-E08; E11-E18",
            "datasets_or_resources": "Human hDRG; human sural nerve; rat DRG; mouse TG; rat resting cornea; mouse wounded cornea",
            "primary_message": "Different targets retain different source components; resting diabetic cornea retains none, whereas wounded diabetic cornea shows provisional residual-component transfer.",
            "source_data": "cross_target_hDRG_component_transportability_tests_2026-08-27.tsv; cross_target_hDRG_component_sample_scores_2026-08-27.tsv",
            "display_files": "cross_target_hDRG_component_transportability_atlas_2026-08-27.png/pdf",
            "evidence_label": "POST-JCI EXPLANATORY ATLAS; RODENT/OCULAR PASSES PROVISIONAL",
            "mandatory_caveat": "Small cohorts and coarse exact-permutation resolution preclude confirmatory cross-species claims.",
            "completion_status": "CORE HEATMAP READY; LEGEND/CONTEXT PANEL PENDING",
        },
        {
            "figure_id": "Figure 6",
            "result_section": "Accessible-compartment boundary",
            "panels": "A, two-cohort PBMC stage test; B, clinical-trial baseline and longitudinal tests; C, tear coverage stop; D, tissue-transport versus biomarker-portability model",
            "evidence_ids": "E19-E21",
            "datasets_or_resources": "GSE95849; GSE185011; GSE302658; PXD062366",
            "primary_message": "Adequately powered blood analyses failed their frozen gates, and tear coverage was insufficient, separating neural-tissue transport from accessible biomarker portability.",
            "source_data": "human_PBMC_stage_projection_meta_analysis_2026-08-27.tsv; human_PBMC_stage_projection_tests_2026-08-27.tsv; GSE302658_clinical_signature_tests_2026-08-27.tsv; PXD062366_hDRG_signature_coverage_2026-08-27.tsv",
            "display_files": "human_PBMC_stage_signature_validation_2026-08-27.png; GSE302658_hDRG_severity_clinical_validation_2026-08-27.png; new tear/model panels required",
            "evidence_label": "PRESPECIFIED FALSIFICATION/BOUNDARY ANALYSES",
            "mandatory_caveat": "A failed projection is not evidence that blood or tears lack all DPN biology; it only rejects portability of the frozen neural programs under the specified gates.",
            "completion_status": "CORE PANELS READY; COMPOSITE PENDING",
        },
    ]
    return pd.DataFrame(rows)


def build_claim_registry() -> tuple[pd.DataFrame, list[str]]:
    signatures = read_tsv(f"hDRG_frozen_primary_stage_signatures_{DATE}.tsv")
    components = read_tsv(f"hDRG_source_defined_transport_components_{DATE}.tsv")
    celltypes = read_tsv(f"hDRG_major_celltype_pairwise_proportion_tests_{DATE}.tsv")
    independent = read_tsv(f"independent_human_DPN_bulk_signature_tests_{DATE}.tsv")
    cross = read_tsv(f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv")
    jci_gate = read_tsv(f"JCI184075_hDRG_component_gate_table_{DATE}.tsv")
    jci_sample = read_tsv(f"JCI184075_hDRG_component_sample_tests_{DATE}.tsv")
    pbmc_meta = read_tsv(f"human_PBMC_stage_projection_meta_analysis_{DATE}.tsv")
    clinical = read_tsv(f"GSE302658_clinical_signature_tests_{DATE}.tsv")
    tear = read_tsv(f"PXD062366_hDRG_signature_coverage_{DATE}.tsv")
    evidence = read_tsv(f"phase0_6_evidence_atlas_{DATE}.tsv")

    assertions: list[str] = []
    signature_counts = signatures.groupby("contrast_id").size().to_dict()
    expected_signatures = {
        "early_allcell_diabetes_vs_control": 656,
        "late_allcell_DPN_vs_diabetes": 3173,
        "late_neuron_DPN_vs_diabetes": 744,
        "severity_neuron_modhigh_vs_low_nageotte": 509,
        "xenium_DPN_vs_control": 60,
    }
    assert signature_counts == expected_signatures
    assertions.append("five frozen source-signature counts match 656/3173/744/509/60")

    component_counts = components.groupby(["module_id", "direction"]).size().unstack(fill_value=0)
    expected_components = {
        "late_shared_concordant_neuronal_core": (77, 119),
        "late_allcell_residual": (1339, 1612),
        "late_neuron_residual": (215, 307),
        "severity_neuron_shared_concordant_core": (104, 33),
        "severity_neuron_residual": (230, 118),
    }
    for component_id, (up, down) in expected_components.items():
        assert int(component_counts.loc[component_id, "up"]) == up
        assert int(component_counts.loc[component_id, "down"]) == down
    assert int((components["module_id"] == "late_directionally_opposed_overlap").sum()) == 26
    assert int((components["module_id"] == "severity_directionally_opposed_overlap").sum()) == 24
    assertions.append("late and severity component counts match the source-only decomposition")

    rows: list[dict[str, str]] = []

    def add(
        claim_id: str,
        result_section: str,
        figures: str,
        evidence_ids: str,
        claim: str,
        source_file: str,
        selector: str,
        extracted: dict[str, object],
        evidence_tier: str,
        caveat: str,
    ) -> None:
        rows.append(
            {
                "claim_id": claim_id,
                "result_section": result_section,
                "figure_panels": figures,
                "evidence_ids": evidence_ids,
                "manuscript_ready_claim": claim,
                "source_file": source_file,
                "row_selector": selector,
                "extracted_values": compact(extracted),
                "evidence_tier": evidence_tier,
                "mandatory_caveat": caveat,
                "audit_status": "PASS",
            }
        )

    add(
        "N001", "Source programs", "Fig. 1B", "E01",
        "The frozen source comprised 656 early all-cell, 3,173 late all-cell, 744 late-neuron, 509 Nageotte-severity and 60 Xenium genes.",
        f"hDRG_frozen_primary_stage_signatures_{DATE}.tsv", "group by contrast_id",
        expected_signatures, "DISCOVERY_ONLY", "The source is a non-peer-reviewed preprint."
    )
    add(
        "N002", "Source programs", "Fig. 1C-D", "E01",
        "Source-only decomposition yielded a 196-gene shared late neuronal core, 522-gene neuron residual and 2,951-gene all-cell residual; 26 directionally opposed overlap genes were audit-only.",
        f"hDRG_source_defined_transport_components_{DATE}.tsv", "group by component_id and direction",
        {"shared_core": 196, "neuron_residual": 522, "allcell_residual": 2951, "opposed_overlap": 26},
        "DISCOVERY_ONLY", "Membership was fixed before target testing."
    )

    cell_rows = []
    for cell_type in ["Neurons", "Fibroblast", "SMC"]:
        row = one(celltypes, cell_type=cell_type, group1="Diabetic", group2="DPN")
        cell_rows.append((cell_type, row))
    add(
        "N003", "Source context", "Fig. 2D", "E22",
        "In source hDRG, neurons decreased from diabetes to DPN (Q=0.0103), whereas fibroblasts (Q=0.0492) and smooth-muscle cells (Q=0.00475) increased.",
        f"hDRG_major_celltype_pairwise_proportion_tests_{DATE}.tsv", "group1=Diabetic; group2=DPN; cell_type in Neurons/Fibroblast/SMC",
        {item[0]: {"q": item[1]["p.adj"], "direction": item[1]["source_article_narrative_for_DPN"]} for item in cell_rows},
        "SOURCE_CONTEXT", "Directions come from the source article narrative and are not independent validation."
    )

    severity = one(independent, contrast_id="severity_neuron_modhigh_vs_low_nageotte")
    assert truth(severity["frozen_gate_supportive"])
    add(
        "N004", "Independent human hDRG", "Fig. 2A-B", "E02",
        "The original Nageotte-severity signature passed in 5 DPN and 7 control hDRGs (gene concordance=0.516, matched Q=0.000400; score difference=0.00809, Hedges g=0.863, sex-stratified exact Q=0.0933), with all leave-one-donor-out differences positive.",
        f"independent_human_DPN_bulk_signature_tests_{DATE}.tsv", "contrast_id=severity_neuron_modhigh_vs_low_nageotte",
        {key: severity[key] for key in ["gene_concordance", "gene_matched_bh_q", "donor_score_difference", "donor_hedges_g", "donor_exact_bh_q", "loo_all_positive"]},
        "INDEPENDENT_HUMAN_VALIDATION", "The cohort lacks diabetes without neuropathy and is age imbalanced."
    )

    human_components = cross[cross["target_id"] == "independent_human_hDRG"].set_index("module_id")
    add(
        "N005", "Independent human hDRG", "Fig. 2C", "E03-E05",
        "In post-primary component analysis, the shared core and all-cell residual were supportive in independent hDRG (g=1.164 and 0.932), whereas the neuron residual was not (g=-0.393).",
        f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv", "target_id=independent_human_hDRG",
        {module: {"g": human_components.loc[module, "hedges_g"], "support": human_components.loc[module, "component_supportive"]} for module in ["late_shared_concordant_neuronal_core", "late_allcell_residual", "late_neuron_residual"]},
        "POST_PRIMARY_EXPLANATORY", "This analysis was performed after the JCI result was known."
    )

    distal = one(jci_gate, target_id="DPN_vs_control", module_id="late_neuron_residual")
    assert truth(distal["component_pass"])
    add(
        "N006", "Human sural nerve", "Fig. 3B-C", "E06",
        "In 6 DPN and 6 control sural nerves, the late neuron residual passed the frozen dual-layer gate (137 mapped up and 252 mapped down genes; gene Q=0.00520; score difference=0.02142, g=1.334, exact Q=0.0281), with all leave-one-out differences positive.",
        f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", "target_id=DPN_vs_control; module_id=late_neuron_residual",
        {key: distal[key] for key in ["mapped_up_n_gene", "mapped_down_n_gene", "gene_concordance", "matched_null_bh_q", "score_mean_difference", "hedges_g", "exact_bh_q", "loo_all_positive"]},
        "FROZEN_INDEPENDENT_HUMAN_TARGET", "Bulk sural nerve is multicellular and does not prove neuronal purity or propagation."
    )

    tpm = one(jci_sample, target_id="DPN_vs_control", expression_source="quantile_normalized_TPM_sensitivity", module_id="late_neuron_residual")
    add(
        "N007", "Human sural nerve", "Fig. 3D", "E06",
        "Quantile-normalized TPM sensitivity analysis reproduced the distal-neuron residual signal (score difference=0.03164, g=1.243, exact Q=0.0476), with all leave-one-out differences positive.",
        f"JCI184075_hDRG_component_sample_tests_{DATE}.tsv", "target_id=DPN_vs_control; expression_source=quantile_normalized_TPM_sensitivity; module_id=late_neuron_residual",
        {key: tpm[key] for key in ["score_mean_difference", "hedges_g", "exact_bh_q", "loo_all_positive"]},
        "SENSITIVITY", "This uses the same donors and is not an independent cohort."
    )

    jci_negative = jci_gate[(jci_gate["target_id"] == "DPN_vs_control") & jci_gate["module_id"].isin(["late_shared_concordant_neuronal_core", "late_allcell_residual"])].set_index("module_id")
    add(
        "N008", "Human sural nerve", "Fig. 3B", "E07-E08",
        "The shared neuronal core failed distally (g=-0.063; gene Q=0.922), and the all-cell residual reversed (g=-1.209; gene Q=1.000).",
        f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", "target_id=DPN_vs_control; two non-supportive late components",
        {module: {"g": jci_negative.loc[module, "hedges_g"], "gene_q": jci_negative.loc[module, "matched_null_bh_q"]} for module in jci_negative.index},
        "FROZEN_NEGATIVE_CONTROLS", "Negative modules demonstrate selectivity rather than absence of all DPN biology."
    )

    axonal = one(jci_gate, target_id="severe_vs_moderate_axonal_loss", module_id="severity_neuron_shared_concordant_core")
    assert truth(axonal["component_pass"])
    add(
        "N009", "Human sural nerve", "Fig. 3E", "E09",
        "Within DPN, the severity-shared neuronal core passed in 13 severe versus 4 moderate axonal-loss specimens (85 mapped up and 22 mapped down genes; gene Q=0.0236; score difference=0.01767, g=1.674, exact Q=0.00504), with all leave-one-out differences positive.",
        f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", "target_id=severe_vs_moderate_axonal_loss; module_id=severity_neuron_shared_concordant_core",
        {key: axonal[key] for key in ["mapped_up_n_gene", "mapped_down_n_gene", "gene_concordance", "matched_null_bh_q", "score_mean_difference", "hedges_g", "exact_bh_q", "loo_all_positive"]},
        "FROZEN_WITHIN_DPN_HUMAN_TARGET", "The moderate group has n=4; the association is cross-sectional."
    )

    severity_residual = one(jci_gate, target_id="severe_vs_moderate_axonal_loss", module_id="severity_neuron_residual")
    assert not truth(severity_residual["component_pass"])
    add(
        "N010", "Human sural nerve", "Fig. 3E", "E10",
        "The severity residual remained a frozen-gate failure despite a positive effect (gene Q=0.0530; g=0.573; exact Q=0.155).",
        f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", "target_id=severe_vs_moderate_axonal_loss; module_id=severity_neuron_residual",
        {key: severity_residual[key] for key in ["matched_null_bh_q", "hedges_g", "exact_bh_q", "component_pass"]},
        "FROZEN_NEGATIVE_NEAR_MISS", "The near-miss cannot be promoted to support."
    )

    for claim_id, target_id, label, evidence_ids, claim, caveat in [
        ("N011", "GSE176017_rat_DRG", "rat DRG", "E11-E13", "In rat painful-DPN progression, the shared core (g=0.797) and all-cell residual (g=1.277) were provisionally supportive, whereas the neuron residual reversed (g=-1.180).", "Only 15 label permutations were possible."),
        ("N012", "GSE329119", "mouse TG", "E14-E15", "In resting diabetic mouse TG, the neuron residual (g=2.273; gene Q=0.0115) and all-cell residual (g=0.820; gene Q=0.000300) were provisionally supportive.", "The cohort was 3 versus 3 and analysis was post-JCI explanatory."),
        ("N013", "GSE227165", "resting rat cornea", "E16", "Resting diabetic rat cornea supported none of the three components (g=-0.531, -0.248 and -0.069 for shared, neuron-residual and all-cell-residual modules).", "Processed libraries, not individual upstream corneas, were the inference units."),
        ("N014", "GSE180490", "wounded mouse cornea", "E17-E18", "In wounded diabetic mouse cornea, the neuron residual (g=0.743; gene Q=0.0444) and all-cell residual (g=0.967; gene Q=0.0255) were provisionally supportive, whereas the shared core failed.", "This is injury-context transfer and not resting-diabetes replication."),
    ]:
        current = cross[cross["target_id"] == target_id].set_index("module_id")
        add(
            claim_id, "Cross-target atlas", "Fig. 5", evidence_ids, claim,
            f"cross_target_hDRG_component_transportability_tests_{DATE}.tsv", f"target_id={target_id}",
            {module: {"g": current.loc[module, "hedges_g"], "gene_q": current.loc[module, "gene_matched_q"], "matched_signature_q": current.loc[module, "matched_signature_q"], "support": current.loc[module, "component_supportive"]} for module in current.index},
            "PROVISIONAL_EXPLANATORY" if target_id != "GSE227165" else "EXPLANATORY_NEGATIVE",
            caveat,
        )

    late_pbmc = one(pbmc_meta, transition="late_DM_to_DPN")
    add(
        "N015", "Accessible compartments", "Fig. 6A", "E19",
        "Across two PBMC cohorts, the frozen late-stage projection was not consistently positive (fixed-effect g=-0.861; I-squared=82.8%), and the PBMC stage gate failed.",
        f"human_PBMC_stage_projection_meta_analysis_{DATE}.tsv", "transition=late_DM_to_DPN",
        {key: late_pbmc[key] for key in ["fixed_effect_hedges_g", "I2_percent", "all_cohort_score_differences_positive", "at_least_one_full_support"]},
        "FROZEN_BOUNDARY", "The two accessions were from the same contributor group and participant non-overlap was not assumed."
    )

    baseline = one(clinical, test_family="primary", contrast_id="severity_neuron_modhigh_vs_low_nageotte", window="baseline", endpoint="npsi_total_score")
    longitudinal = one(clinical, test_family="longitudinal_primary", contrast_id="severity_neuron_modhigh_vs_low_nageotte", window="change", endpoint="npsi_total_score")
    assert not truth(baseline["frozen_gate_pass"]) and not truth(longitudinal["frozen_gate_pass"])
    add(
        "N016", "Accessible compartments", "Fig. 6B", "E20",
        "In 100 baseline PDN blood samples, the severity score was weakly correlated with NPSI (rho=0.070; one-sided permutation P=0.244; matched-signature P=0.280), and the adjusted HC3 exposure test failed (P=0.995).",
        f"GSE302658_clinical_signature_tests_{DATE}.tsv", "test_family=primary; window=baseline; endpoint=npsi_total_score",
        {key: baseline[key] for key in ["n", "rho", "permutation_positive_p", "matched_signature_p", "hc3_exposure_coefficient", "hc3_exposure_p", "frozen_gate_pass"]},
        "FROZEN_CLINICAL_BOUNDARY", "Whole-blood association does not measure ganglion histology."
    )
    add(
        "N017", "Accessible compartments", "Fig. 6B", "E20",
        "Among 95 paired participants, change in the severity score correlated weakly with NPSI change (rho=0.198; one-sided permutation P=0.0273), but the bootstrap lower bound crossed zero and the adjusted Visit-8 exposure test failed (P=0.221); the longitudinal gate therefore failed.",
        f"GSE302658_clinical_signature_tests_{DATE}.tsv", "test_family=longitudinal_primary; window=change; endpoint=npsi_total_score",
        {key: longitudinal[key] for key in ["n", "rho", "permutation_positive_p", "permutation_two_sided_p", "bootstrap_ci_low", "bootstrap_ci_high", "hc3_exposure_coefficient", "hc3_exposure_p", "frozen_gate_pass"]},
        "FROZEN_CLINICAL_BOUNDARY", "Visit 8 occurred under randomized intervention; this is not treatment efficacy or mediation evidence."
    )

    tear_severity = tear[tear["contrast_id"] == "severity_neuron_modhigh_vs_low_nageotte"].set_index("direction")
    add(
        "N018", "Accessible compartments", "Fig. 6C", "E21",
        "The tear proteome retained only 7 severity-up and 33 severity-down proteins at at least 50% detection, failing the prespecified symmetric-coverage requirement before outcome testing.",
        f"PXD062366_hDRG_signature_coverage_{DATE}.tsv", "contrast_id=severity_neuron_modhigh_vs_low_nageotte",
        {direction: {"stable_n": tear_severity.loc[direction, "tear_proteome_gene_n_detection_ge_50pct"], "fraction_original": tear_severity.loc[direction, "fraction_original_detection_ge_50pct"]} for direction in ["up", "down"]},
        "NOT_TESTABLE", "No outcome fishing was performed after the coverage stop."
    )

    assert len(evidence) == 22 and evidence["evidence_id"].is_unique
    assertions.append("evidence atlas contains 22 unique evidence units")
    return pd.DataFrame(rows), assertions


def main() -> None:
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    figure_map = build_figure_map()
    claims, assertions = build_claim_registry()

    figure_path = MANUSCRIPT / f"FIGURE_RESULT_SOURCE_MAP_V1_{DATE}.tsv"
    claim_path = MANUSCRIPT / f"MANUSCRIPT_NUMERICAL_CLAIM_REGISTRY_V1_{DATE}.tsv"
    figure_map.to_csv(figure_path, sep="\t", index=False)
    claims.to_csv(claim_path, sep="\t", index=False)

    source_files = sorted(set(claims["source_file"]))
    audit = {
        "date": DATE,
        "figure_rows": int(len(figure_map)),
        "claim_rows": int(len(claims)),
        "all_claims_pass": bool((claims["audit_status"] == "PASS").all()),
        "assertions": assertions,
        "output_sha256": {
            figure_path.name: sha256(figure_path),
            claim_path.name: sha256(claim_path),
        },
        "source_sha256": {
            name: sha256(TABLES / name)
            for name in source_files
        },
        "boundary": "Registry checks source-table extraction and traceability; it does not convert explanatory or provisional analyses into independent validation.",
    }
    audit_path = MANUSCRIPT / f"MANUSCRIPT_NUMERICAL_CLAIM_AUDIT_V1_{DATE}.json"
    with audit_path.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)

    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
