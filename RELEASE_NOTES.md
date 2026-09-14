# Version 2.0.0

Release date: 2026-09-14

- Replaced legacy symbol handling with corrected canonical GeneID definitions and a strict target-background addendum.
- Retained historical outputs while establishing explicit final replacement precedence.
- Added donor-level spatial ROI, fixed-module deletion, methylation audit and bounded tissue-context analyses.
- Added public clinical extensions for GSE302658, GSE286347, GSE148059/GSE148060/GSE148061, JCI184075 supplements and GSE295206.
- Added complete positive, negative, not-evaluable and blocked evidence registers without pooling P values.
- Updated the manuscript title, citation metadata, environment records and portable path handling.
- Did not redistribute third-party raw data, controlled data, FASTQ, IDAT archives or full public matrices.

This is a major version because corrected program definitions and final source precedence change the numerical reference for the manuscript. The v1.0.0 tag and Zenodo record remain immutable historical artifacts.

## Historical v1.0.0 notes

Initial public reproducibility release for the manuscript:

> Cross-tissue projection of human sensory-ganglion transcriptomic programs in diabetic neuropathy: a multi-cohort computational study

This release provides:

- 21 frozen analysis and audit scripts;
- exact Python environment specifications;
- source-accession, expected-path, and checksum metadata;
- selected processed derivatives required for exact ocular and trigeminal-ganglion projection;
- frozen numerical results, main figures, and traceability records; and
- release-integrity and run-order verification tools.

Third-party raw data are not redistributed. The release reproduces the
source-frozen, non-pooled computational workflow and its reported outputs; it
does not establish causal propagation between tissues or a validated blood or
tear biomarker.

Original software is licensed under MIT. Original non-code content is licensed
under CC BY 4.0 subject to the scope and third-party exclusions in
`LICENSE-CONTENT.md`.
