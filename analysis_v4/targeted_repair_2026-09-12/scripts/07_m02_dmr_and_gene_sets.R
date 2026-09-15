root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root,"05_R_environment/library"); .libPaths(c(lib,.libPaths()))
suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(DMRcate);library(missMethyl)})
outdir <- file.path(root,"07_resource_recovery/M02_results")
anno <- fread(file.path(root,"01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz"),select=c("probe_id","chr","pos"))
members <- fread(file.path(dirname(dirname(root)),"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"))
idcol <- intersect(c("gene_id","GeneID","ncbi_gene_id"),names(members))[1]
pcol <- intersect(c("program","program_id","module_id","Program"),names(members))[1]
stopifnot(!is.na(idcol),!is.na(pcol))
collection <- split(as.character(members[[idcol]]),members[[pcol]])
collection <- lapply(collection,unique)

dmr_status <- list(); gs_status <- list()
for(cohort in c("PROPGER","PROPENG")) {
  e <- fread(file.path(outdir,paste0("M02_",cohort,"_EWAS_ALL_ELIGIBLE.tsv.gz")))
  x <- merge(e,anno,by="probe_id",all.x=TRUE,sort=FALSE)
  x <- x[!is.na(chr)&!is.na(pos)]
  gr <- GRanges(seqnames=x$chr,ranges=IRanges(start=as.integer(x$pos),width=1))
  names(gr) <- x$probe_id
  mcols(gr)$stat <- x$t; mcols(gr)$rawpval <- x$P; mcols(gr)$diff <- x$beta_difference
  mcols(gr)$ind.fdr <- x$q; mcols(gr)$is.sig <- x$q < 0.05
  obj <- new("CpGannotated",ranges=gr,betas=matrix(numeric(0),0,0))
  dmr <- tryCatch(dmrcate(obj,lambda=1000,C=2,pcutoff=0.05,min.cpgs=3),error=function(z) z)
  if(inherits(dmr,"error")) {
    dmr_status[[cohort]] <- data.frame(cohort=cohort,status="NOT_EVALUABLE",reason=conditionMessage(dmr),regions=0)
    writeLines(conditionMessage(dmr),file.path(outdir,paste0("M02_",cohort,"_DMR_ERROR.txt")))
  } else {
    dd <- data.frame(coord=dmr@coord,no_cpgs=dmr@no.cpgs,min_smoothed_fdr=dmr@min_smoothed_fdr,
                     Stouffer=dmr@Stouffer,HMFDR=dmr@HMFDR,Fisher=dmr@Fisher,
                     maxdiff=dmr@maxdiff,meandiff=dmr@meandiff)
    fwrite(dd,file.path(outdir,paste0("M02_",cohort,"_DMR_ALL.tsv")),sep="\t")
    dmr_status[[cohort]] <- data.frame(cohort=cohort,status="COMPLETED_EXPLORATORY",reason="",regions=nrow(dd))
  }
  sig <- e[q<0.05,probe_id]
  gs <- tryCatch(gsameth(sig.cpg=sig,all.cpg=e$probe_id,collection=collection,array.type="EPIC"),error=function(z) z)
  if(inherits(gs,"error")) {
    gs_status[[cohort]] <- data.frame(cohort=cohort,status=if(length(sig)==0) "NOT_EVALUABLE_NO_EWAS_Q05_SEED" else "FAILED_SOFTWARE",reason=conditionMessage(gs),sets=0)
    writeLines(conditionMessage(gs),file.path(outdir,paste0("M02_",cohort,"_GSAMETH_ERROR.txt")))
  } else {
    gg <- as.data.frame(gs); gg$set_id <- rownames(gg); rownames(gg)<-NULL
    fwrite(gg,file.path(outdir,paste0("M02_",cohort,"_GSAMETH_ALL.tsv")),sep="\t")
    gs_status[[cohort]] <- data.frame(cohort=cohort,status="COMPLETED_THRESHOLD_METHOD",reason="",sets=nrow(gg))
  }
}
fwrite(rbindlist(dmr_status),file.path(outdir,"M02_DMR_STATUS.tsv"),sep="\t")
fwrite(rbindlist(gs_status),file.path(outdir,"M02_FIXED_GENE_SET_STATUS.tsv"),sep="\t")
sink(file.path(root,"05_R_environment/M02_sessionInfo.txt"));print(sessionInfo());sink()
