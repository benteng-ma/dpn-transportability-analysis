# Public source inputs

Raw or publisher-supplied third-party files are not redistributed in this archive. Download them from the landing pages in `../../SOURCE_DATA_MANIFEST.tsv` and preserve the exact filenames and relative paths listed there.

Expected top-level input directories are:

```text
data/raw/GSE176017/
data/raw/human_DPN_bulk_PMC8933403/
data/raw/human_hDRG_preprint/
data/raw/human_PBMC_stage_cohorts/GSE95849/
data/raw/human_PBMC_stage_cohorts/GSE185011/
data/raw/human_PDN_trial_GSE302658/
data/raw/human_sural_nerve_JCI184075/supplementary/
data/raw/NCBI_orthology_2026-08-27/
data/raw/PXD062366/
data/raw/reference/
```

GSE176017 per-sample matrices must be extracted into `data/raw/GSE176017/extracted/`. Do not use the quarantined invalid HTML downloads previously returned for JCI supplements; use the valid XLSX files from the official supplementary archive. The reference GMT filenames must match the manifest. Before analysis, run:

```text
python workflow/verify_release.py --check-raw
```

A changed upstream file may legitimately have a different checksum if a repository has replaced or regenerated it. In that case, document the new download date and compare its internal content before treating the analysis as an exact reproduction.
