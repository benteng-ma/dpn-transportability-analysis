# DPN sensory-program transportability analysis

This repository is the versioned computational companion for:

> Anatomical context reshapes sensory ganglion transcriptomic programs in diabetic peripheral neuropathy

Authors: Benteng Ma and Baihua Chen.

## Current manuscript-matched release

Version **2.0.0** corresponds to the corrected Scientific Reports submission candidate dated 2026-09-14. It supersedes v1.0.0 for interpretation of the submitted manuscript but does not delete or rewrite the historical release.

The v2 analysis has three auditable layers:

- `analysis_v3/gene_definition_repair_v1/`: canonical GeneID repair, strict-background addendum, final program membership, coverage, scores and benchmark comparisons.
- `analysis_v4/`: targeted background/module repairs, whole-gene context analyses, methylation method audit, fixed-region follow-up and donor-level spatial ROI closeout.
- `analysis_v5/public_clinical_extensions_2026-09-14/`: public clinical extensions in GSE302658, GSE286347, GSE148059/GSE148060/GSE148061, JCI184075 supplements and GSE295206, plus cross-source evidence integration.

The original v1 workflow remains under `analysis/`, `data/`, `metadata/`, `results/` and `provenance/`. Its immutable release is [v1.0.0](https://github.com/benteng-ma/dpn-transportability-analysis/releases/tag/v1.0.0), archived at [10.5281/zenodo.22151890](https://doi.org/10.5281/zenodo.22151890). That historical DOI must not be cited as the numerical archive for the corrected manuscript.

## Main v2 findings represented here

- Five fixed programs differed between six DPN and six control sural nerves after correction, all opposite to their source direction. Applicable fixed-module deletions did not reverse those effects.
- Seven of nine evaluable programs differed between author-labelled Nageotte and adjacent neuronal regions in six spatial-DRG donors. This is donor-level pathological localization in donors with diabetes history, not independent proof of DPN diagnosis or cell specificity.
- In GSE302658, none of 64 baseline program–symptom tests and none of 64 randomized-treatment interaction tests survived correction. Eight of 64 within-participant symptom-change associations reached BH q<0.05. These are concurrent exploratory associations, not prediction, treatment mediation or causality.
- Expanded blood-methylation, fixed-region and donor-level Nageotte-burden analyses did not provide cross-source corrected support. Not-evaluable analyses remain labelled as such rather than being counted as negative results.

## Reproducibility modes

### Audit the frozen public release

Create the Python environment and run:

```text
python workflow/verify_release_v2.py
python workflow/run_pipeline.py --dry-run
```

The v2 verifier checks required artifacts, Python syntax, manifest integrity, local-path hygiene, file-size limits and the semantic test records supplied with the clinical extension. It does not claim a clean-environment rerun of every upstream raw-data workflow.

### Re-run analyses from public source data

1. Obtain source data from the repositories in `VERSION_2_SOURCE_REGISTER.tsv` and the original `SOURCE_DATA_MANIFEST.tsv`.
2. Set `DPN_PROJECT_ROOT` to the repository root when a script needs inputs outside its immediate analysis directory.
3. Place public inputs in the relative locations documented by the corresponding analysis lock, input audit or script.
4. Use a working copy: some historical execution scripts write same-named derived outputs.
5. Compare newly generated records with the frozen tables and SHA256 manifests.

Large public matrices, FASTQ files, IDAT archives, Visium H5/CLOUPE files and third-party article supplements are not redistributed. Some historical scripts require those public inputs and are retained as provenance rather than represented as one-command portable workflows.

## Important correction history

The v1 source programs contained legacy symbol-to-GeneID ambiguity. The v3.1 repair used standard GeneIDs, preserved genuine cross-source direction conflicts, fixed membership before corrected association reruns, and introduced a strict target-background addendum. Because previous results were already known, this is explicitly a post-result correction and impact analysis—not prospective registration or new blind validation.

The final source authority is:

```text
analysis_v3/gene_definition_repair_v1/strict_background_addendum/consolidated_results/
```

Where a consolidated replacement exists it governs v2. Other historical files remain for audit and must not be mixed with the final strict branch. P9 remains not evaluable for complete RNA scoring because its source-down arm does not meet the fixed coverage rule.

## Public clinical extension

The compact manuscript-matched archive is located at:

```text
analysis_v5/public_clinical_extensions_2026-09-14/
```

It contains analysis-family registers, compact machine-readable results, environment records, evidence matrices, figure-source tables and the scripts used for the public-data extension. No patient crosswalk was guessed, no controlled-access data were used and P values were not pooled across sources.

## Software and environment

- Python dependencies for the original workflow are listed in `requirements.txt` and `environment.yml`.
- The clinical-extension environment snapshot is under `analysis_v5/public_clinical_extensions_2026-09-14/00_admin/environment/`.
- R session information for the corrected and methylation/spatial branches is retained with the corresponding analysis records.
- Public scripts use repository-relative paths by default and accept `DPN_PROJECT_ROOT` where a root override is necessary.

## Licences and data boundaries

Original software is licensed under the MIT License. Original documentation, author-generated metadata, result tables, figures and provenance records are licensed under CC BY 4.0 as scoped in `LICENSE-CONTENT.md`. Third-party source data retain their original licences and are not redistributed or relicensed here.

This repository supports computational audit and bounded reanalysis. It does not establish causal propagation across tissues, clinical prediction, treatment-selection utility, a cell-specific mechanism or independent DPN validation of the spatial reference.

## Citation

Citation metadata are provided in `CITATION.cff`. Cite the version-specific v2.0.0 Zenodo archive at [10.5281/zenodo.22761152](https://doi.org/10.5281/zenodo.22761152). The v1 DOI remains a historical record and should not be substituted for the corrected manuscript-matched release.
