#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: Rscript 04_reproduce_dmr.R CLOSEOUT_DIR")
cdir <- args[1]
p <- dirname(dirname(cdir))
t <- file.path(p,"analysis_v4/targeted_repair_2026-09-12")
.libPaths(c(file.path(t,"05_R_environment/library"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(DMRcate)})
out <- file.path(cdir,"02_M02_audit/DMR_exact_reproduction_v2")
if(file.exists(out)) stop("Refusing to overwrite DMR reproduction directory")
dir.create(out,recursive=TRUE)
anno <- fread(file=file.path(t,"01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz"),select=c("probe_id","chr","pos"))
receipts <- list()
for(cohort in c("PROPGER","PROPENG")) {
  e <- fread(file=file.path(t,"07_resource_recovery/M02_results",paste0("M02_",cohort,"_EWAS_ALL_ELIGIBLE.tsv.gz")))
  x <- merge(e,anno,by="probe_id",all.x=TRUE,sort=FALSE); x <- x[!is.na(chr)&!is.na(pos)]
  gr <- GRanges(seqnames=x$chr,ranges=IRanges(start=as.integer(x$pos),width=1)); names(gr)<-x$probe_id
  mcols(gr)$stat<-x$t; mcols(gr)$rawpval<-x$P; mcols(gr)$diff<-x$beta_difference
  mcols(gr)$ind.fdr<-x$q; mcols(gr)$is.sig<-x$q<0.05
  obj <- new("CpGannotated",ranges=gr,betas=matrix(numeric(0),0,0))
  saveRDS(obj,file.path(out,paste0(cohort,"_CpGannotated.rds")))
  d <- dmrcate(obj,lambda=1000,C=2,pcutoff=0.05,min.cpgs=3)
  got <- data.frame(coord=d@coord,no_cpgs=d@no.cpgs,min_smoothed_fdr=d@min_smoothed_fdr,
                    Stouffer=d@Stouffer,HMFDR=d@HMFDR,Fisher=d@Fisher,
                    maxdiff=d@maxdiff,meandiff=d@meandiff)
  ref <- fread(file=file.path(t,"07_resource_recovery/M02_results",paste0("M02_",cohort,"_DMR_ALL.tsv")))
  write.table(got,file.path(out,paste0(cohort,"_DMR_REPRODUCED.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  exact <- identical(as.character(got$coord),as.character(ref$coord)) && nrow(got)==nrow(ref)
  numeric_max <- max(abs(as.matrix(got[,setdiff(names(got),"coord")])-as.matrix(ref[,setdiff(names(ref),"coord"),with=FALSE])),na.rm=TRUE)
  receipts[[cohort]] <- data.frame(cohort=cohort,n_input=nrow(x),n_q05=sum(x$q<.05),n_regions=nrow(got),coordinates_exact=exact,max_numeric_abs_difference=numeric_max)
}
write.table(do.call(rbind,receipts),file.path(out,"DMR_REPRODUCTION_AUDIT.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(c("Executed parameters: dmrcate(CpGannotated, lambda=1000, C=2, pcutoff=0.05, min.cpgs=3).",
             "The CpGannotated object was reconstructed from the archived complete EWAS tables using t, P, beta difference, BH q and is.sig=(q<0.05).",
             "This audit does not reinterpret absence of single-CpG q<0.05 as proof that region calling is invalid."),file.path(out,"METHOD_NOTE.txt"))
cat("DMR reproduction completed\n")
