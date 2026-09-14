#!/usr/bin/env Rscript
# Capture INSTALLED function bodies/formals; no threshold search or DMR calling.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1L || length(args)>2L) stop("Usage: Rscript audit_methylation_runtime.R NEW_OUT [CpGannotated.rds]")
out <- args[1]
if (file.exists(out)) stop("Output already exists; choose a new audit directory")
dir.create(out, recursive=TRUE)
writeLines(capture.output(sessionInfo()), file.path(out,"SESSION_INFO.txt"))
pkgs <- c("limma","edgeR","minfi","missMethyl","DMRcate","IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
          "IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
versions <- data.frame(package=pkgs, version=NA_character_, available=FALSE)
for (i in seq_along(pkgs)) {
  versions$available[i] <- requireNamespace(pkgs[i], quietly=TRUE)
  if (versions$available[i]) versions$version[i] <- as.character(packageVersion(pkgs[i]))
}
write.table(versions,file.path(out,"PACKAGE_VERSIONS.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
funcs <- list(DMRcate=c("cpg.annotate","dmrcate","extractRanges","changeFDR"),
              missMethyl=c("gsameth","gsaregion"))
for (pkg in names(funcs)) {
  if (!requireNamespace(pkg,quietly=TRUE)) next
  for (fn in funcs[[pkg]]) {
    fun <- get(fn, envir=asNamespace(pkg), inherits=FALSE)
    txt <- c(paste("PACKAGE",pkg,"VERSION",as.character(packageVersion(pkg))),
             "FORMALS",capture.output(dput(formals(fun))),"BODY",deparse(body(fun)))
    writeLines(txt,file.path(out,paste0(pkg,"_",fn,".R.txt")))
  }
}
if (length(args)==2L) {
  if (!requireNamespace("GenomicRanges",quietly=TRUE)) stop("Need GenomicRanges to inspect serialized CpG object")
  obj <- readRDS(args[2])
  writeLines(capture.output(str(obj, max.level=2)),file.path(out,"ANNOTATED_OBJECT_STRUCTURE.txt"))
  if (!isS4(obj) || !"ranges" %in% slotNames(obj)) stop("Unknown CpG object schema: inspect, do not guess")
  rng <- slot(obj,"ranges")
  md <- as.data.frame(S4Vectors::mcols(rng))
  writeLines(names(md),file.path(out,"ANNOTATED_COLUMNS.txt"))
  if ("is.sig" %in% names(md)) {
    write.table(data.frame(n_probes=nrow(md),n_is_sig=sum(md$is.sig %in% TRUE),
                           n_is_sig_missing=sum(is.na(md$is.sig))),
                file.path(out,"ANNOTATED_SIG_COUNTS.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  }
  if ("ind.fdr" %in% names(md)) {
    q <- md$ind.fdr
    if (any(!is.na(q) & (q<0 | q>1))) stop("Invalid stored ind.fdr")
    write.table(data.frame(threshold=c(.05,.1),n_below=c(sum(q<.05,na.rm=TRUE),sum(q<.1,na.rm=TRUE))),
                file.path(out,"ANNOTATED_Q_COUNTS.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  }
}
writeLines("This captures installed defaults and stored object fields, not the actual historical call arguments. Retrieve executed scripts, locks and logs separately.",file.path(out,"IMPORTANT.txt"))
cat("Runtime inspection only; no biological result was recomputed.\n")
