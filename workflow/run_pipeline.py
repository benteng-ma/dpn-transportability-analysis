#!/usr/bin/env python3
"""Print or execute the frozen public-release analysis order."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "analysis" / "scripts"

STEPS = [
    (1, "01_audit_and_map_pxd062366.py", "tear proteome mapping audit"),
    (2, "02_audit_gse176017_matrices.py", "GSE176017 matrix audit"),
    (3, "03_audit_and_extract_hdrg_supplements.py", "source hDRG program freezing"),
    (4, "04_project_hdrg_stages_to_gse176017.py", "rat DRG projection"),
    (5, "05_audit_human_dpn_bulk_supplements.py", "independent hDRG supplement audit"),
    (6, "06_validate_hdrg_signatures_in_independent_human_bulk.py", "independent human hDRG validation"),
    (7, "07_project_hdrg_stages_to_diabetic_tg_cornea.py", "TG and cornea projection"),
    (8, "08_validate_hdrg_stages_in_human_pbmc_cohorts.py", "PBMC boundary analysis"),
    (9, "09_audit_GSE302658_PDN_trial.py", "GSE302658 audit"),
    (10, "10_validate_hdrg_severity_in_GSE302658.py", "clinical whole-blood analysis"),
    (11, "11_audit_PXD062366_hDRG_signature_coverage.py", "tear coverage gate"),
    (12, "12_validate_hdrg_components_in_human_sural_nerve.py", "human sural-nerve validation"),
    (13, "13_build_cross_target_component_transportability_atlas.py", "cross-target atlas"),
    (14, "14_annotate_hdrg_transport_components.py", "functional annotation"),
    (15, "15_extract_hdrg_cell_composition_context.py", "source cell-composition context"),
    (16, "16_build_phase0_6_run_manifest.py", "original-tree run manifest"),
    (17, "17_build_manuscript_traceability_package.py", "claim traceability"),
    (20, "20_build_main_figure_composites.py", "main figure composites"),
    (21, "21_audit_main_figure_composites.py", "main figure audit"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Print the frozen run order by default. Pass --execute only after the "
            "public inputs have been installed and the release has been copied to "
            "a writable working directory."
        )
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--execute", action="store_true", help="Run the selected scripts.")
    action.add_argument("--dry-run", action="store_true", help="Explicitly print without running (default).")
    parser.add_argument("--from-step", type=int, default=1, help="First numbered step to include.")
    parser.add_argument("--to-step", type=int, default=21, help="Last numbered step to include.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = [item for item in STEPS if args.from_step <= item[0] <= args.to_step]
    if not selected:
        raise SystemExit("No runnable public-release steps fall within the requested range.")

    missing_scripts = [name for _, name, _ in selected if not (SCRIPTS / name).is_file()]
    if missing_scripts:
        raise SystemExit("Missing scripts: " + ", ".join(missing_scripts))

    mode = "EXECUTE" if args.execute else "DRY RUN"
    print(f"Release root: {ROOT}")
    print(f"Mode: {mode}")
    if args.execute:
        print("WARNING: same-named derived outputs may be overwritten in this working copy.")

    for number, name, role in selected:
        command = [sys.executable, str(SCRIPTS / name)]
        print(f"[{number:02d}] {role}: {name}", flush=True)
        if args.execute:
            subprocess.run(command, cwd=ROOT, check=True)

    if not args.execute:
        print("No analysis was executed. Add --execute to run the selected steps.")


if __name__ == "__main__":
    main()
