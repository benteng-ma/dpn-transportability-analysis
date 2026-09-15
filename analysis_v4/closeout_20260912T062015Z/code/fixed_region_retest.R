#!/usr/bin/env Rscript
# Retrospective fixed-coordinate follow-up in PROPENG. NOT a new DMR caller.
# Run ONLY after the discovery DMR provenance audit and an explicit local method lock.
# This tests average regional methylation, not DMRcate's omnibus/smoothed hypothesis.
# Inputs: one explicitly constructed RDS list. See SCHEMAS.md for its contract.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: Rscript fixed_region_retest.R JOB.rds NEW_OUTPUT_DIRECTORY")
for (pkg in c("limma", "statmod")) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop(paste("Required local R package missing:", pkg))
}
need <- function(ok, msg) { if (!isTRUE(ok)) stop(msg, call.=FALSE) }
job <- readRDS(args[1])
needed <- c("beta", "metadata", "regions", "membership", "config")
need(all(needed %in% names(job)), "Job missing required objects")
cfg <- job$config
cfg_fields <- c("locked_utc", "role", "discovery_cohort", "test_cohort", "genome_build",
                "coordinates_verified", "discovery_dmr_audit_passed", "qc_verified",
                "phenotype_columns_verified", "covariates", "covariate_types", "factor_levels",
                "branch", "epsilon", "min_cpgs", "min_probe_fraction", "min_residual_df",
                "ebayes_robust", "ebayes_trend", "both_cohort_results_already_known",
                "input_provenance", "annotation_verified")
need(all(cfg_fields %in% names(cfg)), "Incomplete explicit config; no scientific defaults are guessed")
need(is.character(cfg$locked_utc) && nzchar(cfg$locked_utc), "Actual local lock timestamp required")
need(identical(cfg$role, "RETROSPECTIVE_FIXED_REGION_RETEST"), "Invalid role")
need(identical(cfg$discovery_cohort, "PROPGER") && identical(cfg$test_cohort, "PROPENG"), "Cohort roles cannot be swapped")
need(isTRUE(cfg$both_cohort_results_already_known), "Do not label this as a fresh blind holdout")
for (n in c("coordinates_verified", "discovery_dmr_audit_passed", "qc_verified",
            "phenotype_columns_verified", "annotation_verified")) need(isTRUE(cfg[[n]]), paste("Gate not passed:", n))
need(length(cfg$input_provenance) > 0L, "Need actual input and mapping provenance")
need(cfg$branch %in% c("SEX_ADJUSTED_EXPLORATORY", "AGE_BATCH_ADJUSTED"), "Explicit branch required")
need(is.numeric(cfg$epsilon) && length(cfg$epsilon)==1L && cfg$epsilon>0 && cfg$epsilon<0.01, "Invalid epsilon")
need(is.numeric(cfg$min_cpgs) && length(cfg$min_cpgs)==1L && cfg$min_cpgs>=3, "min_cpgs must be >=3")
need(is.numeric(cfg$min_probe_fraction) && cfg$min_probe_fraction>=0 && cfg$min_probe_fraction<=1, "Invalid coverage requirement")
need(cfg$min_residual_df >= 1L, "Invalid residual-df requirement")
need(is.logical(cfg$ebayes_robust) && length(cfg$ebayes_robust)==1L, "Specify robust")
need(is.logical(cfg$ebayes_trend) && length(cfg$ebayes_trend)==1L, "Specify trend")

B <- job$beta
need(is.matrix(B) && is.numeric(B), "beta must be a numeric matrix, CpGs by samples")
need(!is.null(rownames(B)) && !is.null(colnames(B)), "beta dimnames required")
need(!anyDuplicated(rownames(B)) && !anyDuplicated(colnames(B)), "Duplicated probe or sample IDs")
need(all(is.finite(B)) && all(B>=0 & B<=1), "Input must be finite QC-approved beta values; no imputation")
meta <- job$metadata
need(all(c("sample_id", "cohort", "pain_group") %in% names(meta)), "Metadata schema missing")
need(!anyDuplicated(meta$sample_id) && all(!is.na(meta$sample_id)), "Nonunique sample keys")
need(setequal(meta$sample_id, colnames(B)), "Sample keys must match exactly")
meta <- meta[match(colnames(B), meta$sample_id), , drop=FALSE]
need(all(meta$cohort == "PROPENG"), "Test input includes another cohort")
need(all(meta$pain_group %in% c("painless", "painful")), "Do not infer pain classes or include extra controls")
meta$pain_group <- factor(meta$pain_group, levels=c("painless", "painful"))
need(all(table(meta$pain_group)>=2), "Each group needs at least two participants")
need("sex" %in% cfg$covariates, "Preserve declared sex adjustment")
need(all(cfg$covariates %in% names(meta)), "Required covariate unavailable; branch is blocked")
need(all(cfg$covariates %in% names(cfg$covariate_types)), "Missing covariate type")
if (cfg$branch == "AGE_BATCH_ADJUSTED") {
  need(all(c("age", "technical_batch") %in% cfg$covariates), "Main branch needs verified age and batch, not proxies")
}
for (cv in cfg$covariates) {
  tp <- cfg$covariate_types[[cv]]
  need(tp %in% c("numeric", "factor"), "Only explicit numeric or factor covariates")
  if (tp == "numeric") {
    need(is.numeric(meta[[cv]]) && all(is.finite(meta[[cv]])), paste("Non-finite/non-numeric", cv))
  } else {
    lv <- cfg$factor_levels[[cv]]
    need(length(lv)>=2 && all(meta[[cv]] %in% lv), paste("Unverified factor levels", cv))
    meta[[cv]] <- factor(meta[[cv]], levels=lv)
    need(all(table(meta[[cv]])>0), paste("Unused factor level; resolve design without looking at P:", cv))
  }
}
need(!anyNA(meta[, c("pain_group", cfg$covariates), drop=FALSE]), "Do not silently drop patients with missing covariates")
f <- reformulate(c("pain_group", cfg$covariates))
X <- model.matrix(f, meta)
need(qr(X)$rank==ncol(X), "Rank-deficient design; do not drop confounded terms silently")
need(nrow(X)-ncol(X)>=cfg$min_residual_df, "Insufficient residual degrees of freedom")
coef_name <- "pain_grouppainful"
need(coef_name %in% colnames(X), "Missing painful-minus-painless coefficient")

regions <- job$regions
members <- job$membership
need(all(c("region_id", "chr", "start", "end", "discovery_direction", "genome_build") %in% names(regions)), "Region schema missing")
need(nrow(regions)>0L, "No validated discovery regions; stop rather than inventing a list")
need(!anyDuplicated(regions$region_id), "Duplicated regions")
need(all(regions$start>=1 & regions$end>=regions$start), "Coordinates must be 1-based closed")
need(all(regions$start == floor(regions$start) & regions$end == floor(regions$end)), "Noninteger coordinates")
need(all(regions$genome_build==cfg$genome_build), "Genome build mismatch")
need(all(is.na(regions$discovery_direction) | regions$discovery_direction %in% c(-1,0,1)), "Invalid discovery direction")
need(all(c("region_id", "probe_id", "probe_chr", "probe_position", "genome_build") %in% names(members)), "Membership schema missing")
need(!anyDuplicated(members[, c("region_id", "probe_id")]), "Duplicated region/probe pair")
need(all(members$region_id %in% regions$region_id), "Membership of an unselected region")
need(all(members$genome_build == cfg$genome_build), "Probe annotation build mismatch")
rix <- match(members$region_id, regions$region_id)
need(all(members$probe_chr == regions$chr[rix] & members$probe_position >= regions$start[rix] &
         members$probe_position <= regions$end[rix]), "Probe outside fixed discovery coordinates")
# Membership must derive from the frozen discovery region + official QC/annotation.
# No per-probe direction or validation P is used to filter this membership.
clipped_n <- sum(B<=cfg$epsilon | B>=1-cfg$epsilon)
Bc <- B
Bc[Bc < cfg$epsilon] <- cfg$epsilon
Bc[Bc > 1-cfg$epsilon] <- 1-cfg$epsilon
M <- log2(Bc/(1-Bc))
ids <- regions$region_id
rm <- matrix(NA_real_, nrow=length(ids), ncol=ncol(B), dimnames=list(ids, colnames(B)))
rb <- rm
counts <- data.frame(region_id=ids, n_reference_probes=0L, n_test_probes=0L,
                     measured_fraction=NA_real_, status="NOT_EVALUABLE", reason="", stringsAsFactors=FALSE)
for (i in seq_along(ids)) {
  ref <- unique(members$probe_id[members$region_id == ids[i]])
  use <- intersect(ref, rownames(B))
  counts$n_reference_probes[i] <- length(ref)
  counts$n_test_probes[i] <- length(use)
  counts$measured_fraction[i] <- if(length(ref)) length(use)/length(ref) else NA_real_
  if (length(ref)==0L) { counts$reason[i] <- "no_fixed_discovery_probe_membership"; next }
  if (length(use)<cfg$min_cpgs || length(use)/length(ref)<cfg$min_probe_fraction) {
    counts$reason[i] <- "insufficient_QC_qualified_fixed_probes"; next
  }
  rm[i, ] <- colMeans(M[use, , drop=FALSE])
  rb[i, ] <- colMeans(B[use, , drop=FALSE])
  counts$status[i] <- "COMPLETED"
}
keep <- counts$status=="COMPLETED"
ans <- cbind(regions, counts[, setdiff(names(counts), "region_id"), drop=FALSE])
for (n in c("effect_M", "se_M", "ci_low_M", "ci_high_M", "P_two_sided",
            "BH_full_discovery_family", "effect_beta_adjusted", "effect_beta_unadjusted")) ans[[n]] <- NA_real_
ans$direction_matches_discovery <- NA
ans$branch <- cfg$branch
ans$role <- cfg$role
ans$n_selected_regions <- nrow(regions)
ans$n_tested_regions <- sum(keep)
ans$n_participants <- nrow(meta)
if (any(keep)) {
  fit0 <- limma::lmFit(rm[keep, , drop=FALSE], X)
  # Region averages can be constant conditional on the design. Report, do not rescue.
  viable <- is.finite(fit0$sigma) & fit0$sigma>0
  if (any(!viable)) {
    failed_idx <- which(keep)[!viable]
    keep[failed_idx] <- FALSE
    ans$status[failed_idx] <- "NOT_EVALUABLE"
    ans$reason[failed_idx] <- "zero_or_nonfinite_residual_variance"
  }
}
if (any(keep)) {
  fit <- limma::eBayes(limma::lmFit(rm[keep, , drop=FALSE], X),
                      robust=cfg$ebayes_robust, trend=cfg$ebayes_trend)
  j <- match(coef_name, colnames(X))
  b <- fit$coefficients[, j]
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  df <- rep_len(fit$df.total, length(b))
  p <- fit$p.value[, j]
  crit <- qt(0.975, df=df)
  ans$effect_M[keep] <- b
  ans$se_M[keep] <- se
  ans$ci_low_M[keep] <- b-crit*se
  ans$ci_high_M[keep] <- b+crit*se
  ans$P_two_sided[keep] <- p
  # All selected discovery regions remain the validation family, including NE rows.
  ans$BH_full_discovery_family[keep] <- p.adjust(p, method="BH", n=nrow(regions))
  bfit <- limma::lmFit(rb[keep, , drop=FALSE], X)
  ans$effect_beta_adjusted[keep] <- bfit$coefficients[, j]
  ans$effect_beta_unadjusted[keep] <- rowMeans(rb[keep,meta$pain_group=="painful",drop=FALSE])-rowMeans(rb[keep,meta$pain_group=="painless",drop=FALSE])
  # Direction agreement is descriptive and based on beta effects, not a new P test.
  di <- which(keep & !is.na(ans$discovery_direction) & ans$discovery_direction!=0)
  ans$direction_matches_discovery[di] <- sign(ans$effect_beta_adjusted[di])==ans$discovery_direction[di]
}
ans$n_tested_regions <- sum(keep)
out <- args[2]
need(!file.exists(out), "Refusing to overwrite output path")
dir.create(out, recursive=TRUE)
write.table(ans, file.path(out, "FIXED_REGION_RETEST_ALL.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
write.table(data.frame(region_id=rownames(rm), rm, check.names=FALSE), file.path(out,"REGION_MEAN_M.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
write.table(data.frame(region_id=rownames(rb), rb, check.names=FALSE), file.path(out,"REGION_MEAN_BETA.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
write.table(data.frame(sample_id=meta$sample_id, X, check.names=FALSE), file.path(out,"DESIGN_MATRIX.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(cfg, file.path(out,"EXECUTED_CONFIG.rds"))
writeLines(c(capture.output(sessionInfo()), paste("Clipped input beta entries:",clipped_n),
             "Intervals are unadjusted moderated M-scale intervals.",
             "No new regions called in the test cohort. Known-results retrospective evaluation.",
             "No claim of independent clinical validation from this helper alone."), file.path(out,"RUNTIME_AND_LIMITS.txt"))
cat("Finished fixed-region follow-up. Output:", normalizePath(out), "\n")
