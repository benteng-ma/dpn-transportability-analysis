root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root, "05_R_environment/library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", lib=lib, quiet=FALSE)
}
needed <- c("limma", "edgeR")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  BiocManager::install(missing, lib=lib, ask=FALSE, update=FALSE)
}
stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly=TRUE)))

sink(file.path(root, "05_R_environment/sessionInfo.txt"))
cat("R executable:", file.path(R.home("bin"), "Rscript.exe"), "\n")
cat("Project library:", lib, "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")
cat("edgeR:", as.character(packageVersion("edgeR")), "\n\n")
print(sessionInfo())
sink()

set.seed(20260912)
toy <- matrix(rnbinom(200, mu=50, size=5), nrow=20)
group <- factor(rep(c("A", "B"), each=5))
y <- edgeR::DGEList(toy, group=group)
keep <- edgeR::filterByExpr(y, group=group)
y <- edgeR::calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
design <- model.matrix(~group)
v <- limma::voom(y, design, plot=FALSE)
fit <- limma::eBayes(limma::lmFit(v, design), robust=TRUE)
smoke <- limma::topTable(fit, coef="groupB", number=Inf, sort.by="none")
stopifnot(nrow(smoke) > 0, all(is.finite(smoke$P.Value)))
write.table(data.frame(test="edgeR_limma_voom_smoke", status="PASS", retained_genes=nrow(smoke)),
            file=file.path(root, "05_R_environment/R_SMOKE_TEST.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("R environment and smoke test completed\n")

