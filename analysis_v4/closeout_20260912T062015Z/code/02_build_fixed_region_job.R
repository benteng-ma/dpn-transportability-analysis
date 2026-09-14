#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript 02_build_fixed_region_job.R CLOSEOUT_DIR OUT_RDS")
cdir <- args[1]
read_tsv <- function(x) read.delim(x, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("NA",""))
b <- read_tsv(file.path(cdir,"02_M02_audit/PROPENG_FIXED_REGION_BETA.tsv.gz"))
rownames(b) <- b$probe_id; b$probe_id <- NULL; B <- as.matrix(b); storage.mode(B) <- "double"
meta <- read_tsv(file.path(cdir,"02_M02_audit/PROPENG_FIXED_REGION_METADATA.tsv"))
regions <- read_tsv(file.path(cdir,"02_M02_audit/PROPGER_FIXED_REGIONS.tsv"))
membership <- read_tsv(file.path(cdir,"02_M02_audit/PROPGER_FIXED_REGION_MEMBERSHIP.tsv"))
sha <- function(path) unname(tools::md5sum(path))
cfg <- list(
  locked_utc=format(Sys.time(),tz="UTC",usetz=TRUE), role="RETROSPECTIVE_FIXED_REGION_RETEST",
  discovery_cohort="PROPGER", test_cohort="PROPENG", genome_build="hg19",
  coordinates_verified=TRUE, discovery_dmr_audit_passed=TRUE, qc_verified=TRUE,
  phenotype_columns_verified=TRUE, annotation_verified=TRUE,
  both_cohort_results_already_known=TRUE, covariates=c("sex"),
  covariate_types=c(sex="factor"), factor_levels=list(sex=c("female","male")),
  branch="SEX_ADJUSTED_EXPLORATORY", epsilon=1e-6, min_cpgs=3,
  min_probe_fraction=0, min_residual_df=1, ebayes_robust=TRUE, ebayes_trend=FALSE,
  input_provenance=list(
    beta="official GSE286347 beta matrix; SHA256 recorded in M02_INPUT_PROVENANCE.json",
    regions="all 167 valid PROPGER DMRcate regions from the executed targeted-repair branch",
    membership="PROPGER QC-qualified EPIC probes within fixed hg19 1-based closed coordinates",
    note="Both cohort DMR results were known before this retrospective fixed-region retest"
  )
)
saveRDS(list(beta=B, metadata=meta, regions=regions, membership=membership, config=cfg), args[2])
cat("Saved",args[2],"with",nrow(B),"probes and",ncol(B),"participants\n")
