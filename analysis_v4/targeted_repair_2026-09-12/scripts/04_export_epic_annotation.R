root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root, "05_R_environment/library")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})
a <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
a$probe_id <- rownames(a)
keep <- intersect(c("probe_id","chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Group"), colnames(a))
out <- a[,keep,drop=FALSE]
write.table(out, gzfile(file.path(root,"01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz")), sep="\t", row.names=FALSE, quote=FALSE, na="")
write.table(data.frame(annotation_package="IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                       package_version=as.character(packageVersion("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")),
                       rows=nrow(out), columns=paste(colnames(out),collapse=";")),
            file.path(root,"01_inputs/EPIC_V1_HG19_ANNOTATION_AUDIT.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
