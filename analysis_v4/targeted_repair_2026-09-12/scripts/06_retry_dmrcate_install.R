root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root,"05_R_environment/library")
.libPaths(c(lib,.libPaths())); options(repos=c(CRAN="https://cloud.r-project.org"))
stopifnot(requireNamespace("BiocManager",quietly=TRUE))
ordered <- c("TxDb.Hsapiens.UCSC.hg19.knownGene","FDb.InfiniumMethylation.hg19","DMRcate")
for (p in ordered) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,lib=lib,ask=FALSE,update=FALSE)
s <- data.frame(package=c("minfi","missMethyl","IlluminaHumanMethylationEPICanno.ilm10b4.hg19",ordered),
                available=vapply(c("minfi","missMethyl","IlluminaHumanMethylationEPICanno.ilm10b4.hg19",ordered),requireNamespace,logical(1),quietly=TRUE))
write.table(s,file.path(root,"05_R_environment/METHYLATION_PACKAGES_AFTER_RETRY.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(all(s$available))
