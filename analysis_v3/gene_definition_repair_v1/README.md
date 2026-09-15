# DPN canonical gene-definition repair v1

This is a local, retrospective correction and impact-assessment bundle after known v2 and Batch01-04 results. It is not a submission-ready manuscript, preregistration, independent external validation, or a public release. Original results are preserved.

## Final strict-background addendum

Handoff checks identified four mapping failures in the first target-background implementation: ambiguous official symbols HBD, MMD2 and TEC were assigned by a legacy last-row dictionary. Source members were already ambiguity-safe and unchanged. A separate strict branch was locked before its own outcomes and executed; initial results and failed checks remain preserved. **Use `strict_background_addendum/consolidated_results/` in preference to the same-named `results/` file whenever present.** The final PDF, `MAIN_EFFECT_COVERAGE_COMPARISON.tsv`, and `DONOR_SCORE_COMPARISON.tsv.gz` read that consolidated view. Initial presentation is retained in `initial_presentation_before_strict_addendum/` and is not the final report.

See `reports/STRICT_ADDENDUM_SUMMARY.json`, addendum input lock and implementation-impact tables. The unchanged GSE24290 fits were reused. Final benchmark Holm families retain all 40 eligible tests/background. Historical P changes in this addendum affect one Xenium-parent case-state test only (q remains 1); substantive conclusions are unchanged. This is a disclosed implementation correction, not a new opportunity for favorable threshold or membership choice.

## Read first

- `reports/REPAIR_REVIEW_REPORT.pdf`: merged review, scientific boundaries, figures and primary numeric comparisons.
- `reports/REPAIR_REPORT_CN.md`: Chinese conclusions and limitations.
- `config/REPAIR_RULES.md`, `config/LOCK.json`, `config/LOCKED_FILES.tsv`: source-only rules and actual pre-corrected-association lock.
- `STATUS.tsv`: dependency-by-dependency completion and exclusions.
- `results/MAIN_EFFECT_COVERAGE_COMPARISON.tsv`, `results/DONOR_SCORE_COMPARISON.tsv.gz`, `results/OLD_NEW_MEMBER_ROW_COMPARISON.tsv`: explicit old/new comparisons.
- `results/HISTORICAL_DONOR_TEST_IMPACT.tsv`: historical one-sided tests and original families. Do not substitute common-method P/q for these historical values.
- `results/REAL_BENCHMARK_ALL.tsv`, `results/MATCHED_RANDOM_SPLITS.tsv.gz`: all empirical benchmark tasks and actual random membership draws, not selected favorable results.
- `tests/`, `logs/`, and `MANIFEST_SHA256.tsv`: numeric/integrity checks, actual execution including failures, and payload hashes.

## Interpretation

The severity core has 16 up and 7 down GeneIDs. It fails the unchanged minimum of ten measured genes in each arm. Its old axonal-loss association cannot be carried into corrected claims; this is not evidence of biological absence. Late-program historical associations are not all erased, but the corrected direct/random decomposition benchmark does not establish a multiple-testing-corrected advantage. CIAP reference findings are not DPN validation. Two ALC products are not independent cohorts.

## Reproduction and portability

The local scripts live beneath the original DPN project in `analysis_v3/gene_definition_repair_v1`. They explicitly reuse archived, hashed numerical helper code and source files from that project. The bundle includes new scripts, outputs, matrices and source-row audit, but is NOT a standalone public reanalysis distribution: original workbooks, large cell-reference caches and complete environments remain required in the authorized DPN project. `inputs/SOURCE_FILES.tsv`, `inputs/TARGET_SOURCE_FILES.tsv`, `inputs/PROTECTED_FILES.tsv` and the script dependency inventory identify original inputs. Do not run legacy helper main functions directly against the original project.

Use the existing Batch04 project Python environment, not a global installation. Source extraction/locking scripts 00 and 01 refuse to overwrite an existing lock. Do not rerun them in this frozen directory. For a new reanalysis, create another explicitly named directory and retain this version. Later scripts depend on earlier outputs; their actual execution order is in logs (not an assertion of pre-result registration).

R exclusion sensitivity used R 4.6.1 and the read-only Batch02 R library. A missing executable attempt, an incompatible R/rlang attempt, and a parent-context join-key failure are retained in logs; subsequent successful runs are separately recorded. Successful scientific outputs were not replaced merely to hide failed packaging/join attempts.

Legacy-style files ending in `2026-08-27` retain interface-compatible names only. They were recomputed in this correction run; actual timestamps and source-code hashes are in logs. Reports and result paths are project-relative; private logs may retain runtime paths. This local review bundle is not approved for public redistribution.

PNG/TIFF are the review figures. Machine-readable TSV/GZ files contain full results; there is no Excel-only supplement. A PDF reader table is not a replacement for those full data. No GitHub, Zenodo, DOI, manuscript, author declaration or funding record was changed.
