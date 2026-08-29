# DPN transportability analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22151890.svg)](https://doi.org/10.5281/zenodo.22151890)

This repository is the reproducibility companion for:

> Cross-tissue projection of human sensory-ganglion transcriptomic programs in diabetic neuropathy: a multi-cohort computational study

Authors: Benteng Ma and Baihua Chen.

## Release status

Version **1.0.0** was released on 2026-08-29 and is publicly archived in Zenodo at [https://doi.org/10.5281/zenodo.22151890](https://doi.org/10.5281/zenodo.22151890). The corresponding GitHub release is [v1.0.0](https://github.com/benteng-ma/dpn-transportability-analysis/releases/tag/v1.0.0). For exact reproducibility of the submitted manuscript, cite the version-specific Zenodo DOI rather than the all-versions concept DOI.

## What this release contains

- `analysis/scripts/`: 21 frozen Python analysis and audit scripts.
- `data/processed/ocular/`: six frozen processed derivatives used for the ocular and trigeminal-ganglion projection, plus their checksums and provenance.
- `data/raw/`: an intentionally empty input area with acquisition instructions. Third-party raw data are not redistributed.
- `metadata/`: frozen sample and cohort metadata created during the study.
- `results/tables/`: frozen numerical results and analysis-level QC records.
- `results/figures/`: frozen main figures and the two source figure pairs needed to rebuild the composite set.
- `provenance/`: numerical-claim, bibliographic, dataset, manuscript, and figure audit records.
- `workflow/`: a run-order helper and release validator.

The scientific design is source-frozen and non-pooled: source-program membership and direction were set before target scoring, and every target cohort was processed independently. Matrices and study-level effects were not pooled across datasets.

## Two reproducibility modes

### 1. Audit the frozen release without downloading raw data

Create the environment, then run:

```text
python workflow/verify_release.py
python workflow/run_pipeline.py --dry-run
```

This checks the archive structure, code syntax, local-path hygiene, checksums, frozen results, and expected input layout. It does not repeat upstream analyses that require third-party source files.

### 2. Repeat the full analysis from public source data

1. Recreate the Python environment using `environment.yml` or `requirements.txt`.
2. Obtain the public inputs listed in `SOURCE_DATA_MANIFEST.tsv` from the cited repositories or article supplements.
3. Place them under `data/raw/` using the exact relative paths in the manifest and `data/raw/README.md`.
4. Run `python workflow/verify_release.py --check-raw`.
5. Work in a copy of this release, because the scripts overwrite same-named derived outputs.
6. Run `python workflow/run_pipeline.py --execute`.
7. Run `python workflow/verify_release.py --check-raw` again and compare the generated output hashes or numerical tables with the frozen release.

Some upstream public resources are large. They are deliberately referenced rather than bundled. Repository landing pages, expected paths, source roles, and hashes of the locally analysed copies are recorded in `SOURCE_DATA_MANIFEST.tsv`.

## Frozen execution order

| Step | Script | Role |
|---:|---|---|
| 01 | `01_audit_and_map_pxd062366.py` | Audit and map the tear-proteome source. |
| 02 | `02_audit_gse176017_matrices.py` | Audit GSE176017 animal matrices and create pseudobulk counts. |
| 03 | `03_audit_and_extract_hdrg_supplements.py` | Extract and freeze source hDRG programs. |
| 04 | `04_project_hdrg_stages_to_gse176017.py` | Project source programs into rat DRG. |
| 05 | `05_audit_human_dpn_bulk_supplements.py` | Audit the independent human hDRG supplements. |
| 06 | `06_validate_hdrg_signatures_in_independent_human_bulk.py` | Perform independent human hDRG validation. |
| 07 | `07_project_hdrg_stages_to_diabetic_tg_cornea.py` | Project source programs into TG and corneal datasets. |
| 08 | `08_validate_hdrg_stages_in_human_pbmc_cohorts.py` | Test two PBMC cohorts. |
| 09 | `09_audit_GSE302658_PDN_trial.py` | Audit GSE302658 and construct sample metadata. |
| 10 | `10_validate_hdrg_severity_in_GSE302658.py` | Test the clinical whole-blood cohort. |
| 11 | `11_audit_PXD062366_hDRG_signature_coverage.py` | Apply the prespecified tear coverage gate. |
| 12 | `12_validate_hdrg_components_in_human_sural_nerve.py` | Test components in human sural nerve. |
| 13 | `13_build_cross_target_component_transportability_atlas.py` | Build the cross-target component atlas. |
| 14 | `14_annotate_hdrg_transport_components.py` | Perform competitive functional annotation. |
| 15 | `15_extract_hdrg_cell_composition_context.py` | Extract source cell-composition context. |
| 16 | `16_build_phase0_6_run_manifest.py` | Build the original project manifest. |
| 17 | `17_build_manuscript_traceability_package.py` | Build result-to-claim traceability records. |
| 18 | `18_audit_manuscript_numerical_consistency.py` | Manuscript-development audit; retained for provenance. |
| 19 | `19_audit_full_manuscript.py` | Manuscript-development audit; retained for provenance. |
| 20 | `20_build_main_figure_composites.py` | Rebuild the six main figures. |
| 21 | `21_audit_main_figure_composites.py` | Audit figure dimensions and hashes. |

Steps 18 and 19 refer to the editable manuscript-development files and are not invoked by the public pipeline helper. Their frozen audit outputs are supplied under `provenance/`. Step 16 describes the original full project tree; `workflow/verify_release.py` is the authoritative validator for this public release layout.

## Interpretation boundaries

This archive supports computational reproducibility, not causal inference. Association or signature transport across tissues does not establish propagation from ganglion soma to peripheral axon. Post-primary animal and ocular analyses are explanatory and cannot substitute for independent human validation. Failed blood or tear transfer rejects a direct proxy under the frozen test; it does not establish absence of all diabetes-related biology in those compartments.

## Environment

The release was tested with Python 3.13.11. Exact package versions are recorded in `requirements.txt`; a portable environment recipe is supplied in `environment.yml`. No internet access is used by the analysis scripts once the required source files are present.

## License and citation

Original analysis software is licensed under the MIT License in `LICENSE`. Original documentation, author-generated metadata, frozen result tables, figures, and provenance records are licensed under CC BY 4.0 as scoped in `LICENSE-CONTENT.md`. Third-party source data remain governed by their original repository or publisher terms and are not redistributed or relicensed here. Citation metadata are provided in `CITATION.cff`; the archived version should be cited as DOI [10.5281/zenodo.22151890](https://doi.org/10.5281/zenodo.22151890).

## Contact

Questions about the scientific analysis should be directed to the corresponding author through the contact information in the associated manuscript. No support commitment is implied by this archive.
