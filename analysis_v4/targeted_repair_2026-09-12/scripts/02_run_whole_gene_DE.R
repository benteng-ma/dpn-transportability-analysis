root <- "${DPN_PROJECT_ROOT}/analysis_v4/targeted_repair_2026-09-12"
lib <- file.path(root, "05_R_environment/library")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(limma)})
set.seed(20260912)

read_counts <- function(path) {
  x <- read.delim(path, row.names=1, check.names=FALSE)
  x <- as.matrix(x); storage.mode(x) <- "numeric"
  stopifnot(all(is.finite(x)), all(x >= 0), all(abs(x-round(x)) < 1e-8))
  x
}

write_full <- function(tt, path, gene_ids) {
  stopifnot(setequal(rownames(tt), gene_ids))
  tt$gene_id <- rownames(tt)
  tt <- tt[,c("gene_id", setdiff(colnames(tt), "gene_id"))]
  write.table(tt, path, sep="\t", row.names=FALSE, quote=FALSE, na="")
}

# M04 primary and sex-adjusted sensitivity. Group is inseparable from anatomy and surgical source.
m4 <- read_counts(file.path(root, "01_inputs/M04_DE_COUNTS_SUM.tsv.gz"))
meta4 <- read.delim(file.path(root, "01_inputs/M04_SAMPLE_METADATA.tsv"), check.names=FALSE)
meta4 <- meta4[match(colnames(m4), meta4$sample_id),]
stopifnot(identical(colnames(m4), meta4$sample_id))
meta4$case <- factor(ifelse(meta4$group == "Morton’s neuroma", "Morton", "Control"), levels=c("Control","Morton"))
meta4$sex <- factor(meta4$sex)
y4 <- DGEList(m4)
design4 <- model.matrix(~case, meta4)
keep4 <- filterByExpr(y4, design4)
y4 <- calcNormFactors(y4[keep4,,keep.lib.sizes=FALSE], method="TMM")
v4 <- voom(y4, design4, plot=FALSE)
fit4 <- eBayes(lmFit(v4, design4), robust=TRUE)
tt4 <- topTable(fit4, coef="caseMorton", number=Inf, sort.by="none")
write_full(tt4, file.path(root, "06_whole_gene_DE/results/M04_PRIMARY_DE_ALL_GENES.tsv"), rownames(m4)[keep4])

design4s <- model.matrix(~sex + case, meta4)
keep4s <- filterByExpr(DGEList(m4), design4s)
y4s <- calcNormFactors(DGEList(m4[keep4s,]), method="TMM")
v4s <- voom(y4s, design4s, plot=FALSE)
fit4s <- eBayes(lmFit(v4s, design4s), robust=TRUE)
tt4s <- topTable(fit4s, coef="caseMorton", number=Inf, sort.by="none")
write_full(tt4s, file.path(root, "06_whole_gene_DE/results/M04_SEX_ADJUSTED_DE_ALL_GENES.tsv"), rownames(m4)[keep4s])
write.table(data.frame(module="M04", model=c("primary_group_only","sex_adjusted_sensitivity"), tested_genes=c(nrow(tt4),nrow(tt4s)), q_lt_005=c(sum(tt4$adj.P.Val<0.05),sum(tt4s$adj.P.Val<0.05)), q_lt_010=c(sum(tt4$adj.P.Val<0.10),sum(tt4s$adj.P.Val<0.10))), file=file.path(root,"06_whole_gene_DE/results/M04_DE_SUMMARY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# M05 cell-means model with donor correlation. Disease remains a donor-level observational contrast.
m5 <- read_counts(file.path(root, "01_inputs/M05_DE_COUNTS_SUM.tsv.gz"))
meta5 <- read.delim(file.path(root, "01_inputs/M05_SAMPLE_METADATA.tsv"), check.names=FALSE)
meta5 <- meta5[match(colnames(m5), meta5$sample_id),]
stopifnot(identical(colnames(m5), meta5$sample_id), length(unique(meta5$donor)) == 9)
meta5$cell <- factor(paste(meta5$muscle, meta5$group, sep="_"), levels=c("AH_MNC","AH_DPN","MG_MNC","MG_DPN"))
design5 <- model.matrix(~0+cell, meta5); colnames(design5) <- levels(meta5$cell)
y5 <- DGEList(m5)
keep5 <- filterByExpr(y5, design5)
y5 <- calcNormFactors(y5[keep5,,keep.lib.sizes=FALSE], method="TMM")
v5a <- voom(y5, design5, plot=FALSE)
corfit <- duplicateCorrelation(v5a, design5, block=meta5$donor)
v5 <- voom(y5, design5, plot=FALSE, block=meta5$donor, correlation=corfit$consensus.correlation)
fit5 <- lmFit(v5, design5, block=meta5$donor, correlation=corfit$consensus.correlation)
contr <- makeContrasts(AH_DPN_minus_MNC=AH_DPN-AH_MNC, MG_DPN_minus_MNC=MG_DPN-MG_MNC, diagnosis_by_muscle=(AH_DPN-AH_MNC)-(MG_DPN-MG_MNC), levels=design5)
fit5c <- eBayes(contrasts.fit(fit5, contr), robust=TRUE)
sum5 <- list()
for (coef in colnames(contr)) {
  tt <- topTable(fit5c, coef=coef, number=Inf, sort.by="none")
  out <- file.path(root, "06_whole_gene_DE/results", paste0("M05_",coef,"_DE_ALL_GENES.tsv"))
  write_full(tt, out, rownames(m5)[keep5])
  sum5[[coef]] <- data.frame(module="M05", contrast=coef, tested_genes=nrow(tt), q_lt_005=sum(tt$adj.P.Val<0.05), q_lt_010=sum(tt$adj.P.Val<0.10), duplicate_correlation=corfit$consensus.correlation)
}
write.table(do.call(rbind,sum5), file=file.path(root,"06_whole_gene_DE/results/M05_DE_SUMMARY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

png(file.path(root,"06_whole_gene_DE/figures/M04_M05_DE_diagnostic.png"), width=2400, height=1600, res=300)
par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(tt4$logFC, -log10(pmax(tt4$P.Value, .Machine$double.xmin)), pch=16, cex=.35, col="#4d4d4d", xlab="M04 log2 fold change", ylab="-log10(P)"); abline(h=-log10(.05), lty=2, col="#777777")
tt5 <- topTable(fit5c, coef="AH_DPN_minus_MNC", number=Inf, sort.by="none")
plot(tt5$logFC, -log10(pmax(tt5$P.Value, .Machine$double.xmin)), pch=16, cex=.35, col="#4d4d4d", xlab="M05 AH log2 fold change", ylab="-log10(P)"); abline(h=-log10(.05), lty=2, col="#777777")
dev.off()

sink(file.path(root, "05_R_environment/DE_sessionInfo.txt")); print(sessionInfo()); sink()
cat("Whole-gene DE completed\n")
