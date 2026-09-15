root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root, "05_R_environment/library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", lib=lib)
pkgs <- c("missMethyl", "DMRcate", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) BiocManager::install(missing, lib=lib, ask=FALSE, update=FALSE)
status <- data.frame(package=pkgs,
                     available=vapply(pkgs, requireNamespace, logical(1), quietly=TRUE),
                     version=vapply(pkgs, function(x) if(requireNamespace(x,quietly=TRUE)) as.character(packageVersion(x)) else NA_character_, character(1)))
write.table(status, file.path(root,"05_R_environment/METHYLATION_PACKAGES.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="")
stopifnot(all(status$available))
sink(file.path(root,"05_R_environment/METHYLATION_sessionInfo.txt")); print(sessionInfo()); sink()
