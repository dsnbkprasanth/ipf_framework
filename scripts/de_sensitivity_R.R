#!/usr/bin/env Rscript
# Optional gold-standard DE sensitivity in R: limma-voom + edgeR (RNA-seq, GSE150910)
# and limma (microarray, GSE47460). Requires R with Bioconductor.
#
# Install once (in R):
#   install.packages("BiocManager")
#   BiocManager::install(c("limma","edgeR","DESeq2"))
#
# Run:
#   Rscript scripts/de_sensitivity_R.R
#
# Compares the limma-voom replicated RNA-seq set with the pipeline's primary
# signature (results_real/tables/S_stage1_replicated_signature.csv).

suppressMessages({library(limma); library(edgeR)})

counts_file <- "data/GSE150910/GSE150910_gene-level_count_file.csv.gz"
primary_file <- "results_real/tables/S_stage1_replicated_signature.csv"

cts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
cts <- rowsum(cts, group = rownames(cts))                 # collapse duplicate symbols
grp <- ifelse(grepl("^ipf", tolower(colnames(cts))), "IPF",
        ifelse(grepl("^control", tolower(colnames(cts))), "control", NA))
keep_s <- !is.na(grp); cts <- cts[, keep_s]; grp <- factor(grp[keep_s], levels = c("control","IPF"))
cat(sprintf("RNA-seq: %d genes x %d samples (%d IPF / %d control)\n",
            nrow(cts), ncol(cts), sum(grp=="IPF"), sum(grp=="control")))

# limma-voom
dge <- DGEList(counts = round(as.matrix(cts)))
keep <- filterByExpr(dge, group = grp); dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ grp)
v <- voom(dge, design)
fit <- eBayes(lmFit(v, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
voom_up <- rownames(tt)[tt$adj.P.Val < 0.05 & tt$logFC > 0]
voom_dn <- rownames(tt)[tt$adj.P.Val < 0.05 & tt$logFC < 0]
cat(sprintf("limma-voom (GSE150910) significant: up=%d down=%d\n", length(voom_up), length(voom_dn)))

# overlap with primary signature
prim <- read.csv(primary_file, row.names = 1, check.names = FALSE)
prim_up <- rownames(prim)[prim$replicated_up == "True" | prim$replicated_up == TRUE]
prim_dn <- rownames(prim)[prim$replicated_down == "True" | prim$replicated_down == TRUE]
fdr47 <- grep("^fdr_.*47460", colnames(prim), value = TRUE)[1]
lfc47 <- grep("^log2FC_.*47460", colnames(prim), value = TRUE)[1]
common <- intersect(rownames(tt), rownames(prim))
a_up <- common[prim[common, fdr47] < 0.05 & prim[common, lfc47] > 0]
a_dn <- common[prim[common, fdr47] < 0.05 & prim[common, lfc47] < 0]
hy_up <- intersect(intersect(voom_up, common), a_up)
hy_dn <- intersect(intersect(voom_dn, common), a_dn)
cat(sprintf("limma-voom-hybrid replicated: up=%d down=%d\n", length(hy_up), length(hy_dn)))
cat(sprintf("  %% within primary: up=%.1f%%  down=%.1f%%\n",
            100*length(intersect(hy_up, prim_up))/max(length(hy_up),1),
            100*length(intersect(hy_dn, prim_dn))/max(length(hy_dn),1)))

out <- data.frame(gene = rownames(tt), logFC = tt$logFC, adj_P = tt$adj.P.Val)
write.csv(out, "results_real/tables/S2b_limmavoom_GSE150910.csv", row.names = FALSE)
cat("Wrote results_real/tables/S2b_limmavoom_GSE150910.csv\n")
