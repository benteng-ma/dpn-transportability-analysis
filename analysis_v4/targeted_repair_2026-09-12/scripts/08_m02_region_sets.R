root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root,"05_R_environment/library"); .libPaths(c(lib,.libPaths()))
suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(missMethyl)})
outdir <- file.path(root,"07_resource_recovery/M02_results")
members <- fread(file.path(dirname(dirname(root)),"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"))
collection <- lapply(split(as.character(members$gene_id),members$module_id),unique)
status <- list()
for(cohort in c("PROPGER","PROPENG")) {
  d <- fread(file.path(outdir,paste0("M02_",cohort,"_DMR_ALL.tsv")))
  z <- tstrsplit(d$coord,"[:-]")
  gr <- GRanges(seqnames=z[[1]],ranges=IRanges(as.integer(z[[2]]),as.integer(z[[3]])))
  e <- fread(file.path(outdir,paste0("M02_",cohort,"_EWAS_ALL_ELIGIBLE.tsv.gz")),select="probe_id")
  gs <- tryCatch(gsaregion(regions=gr,all.cpg=e$probe_id,collection=collection,array.type="EPIC"),error=function(x)x)
  if(inherits(gs,"error")) {
    status[[cohort]] <- data.frame(cohort=cohort,status="FAILED_SOFTWARE",sets=0,reason=conditionMessage(gs))
    writeLines(conditionMessage(gs),file.path(outdir,paste0("M02_",cohort,"_GSAREGION_ERROR.txt")))
  } else {
    gg <- as.data.frame(gs); gg$set_id <- rownames(gg); rownames(gg)<-NULL
    pnames <- grep("^P\\.",names(gg),value=TRUE)
    for(pc in pnames) gg[[paste0(pc,"_BH")]] <- p.adjust(gg[[pc]],method="BH")
    fwrite(gg,file.path(outdir,paste0("M02_",cohort,"_GSAREGION_ALL.tsv")),sep="\t")
    status[[cohort]] <- data.frame(cohort=cohort,status="COMPLETED_EXPLORATORY_REGION_SET",sets=nrow(gg),reason="")
  }
}
fwrite(rbindlist(status),file.path(outdir,"M02_REGION_SET_STATUS.tsv"),sep="\t")
