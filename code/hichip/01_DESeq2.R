library(diffloop)
library(matrixStats)
library(dplyr)

cpm <- function(mat){
  t(t(mat)/colSums(mat))*1000000
}

# Collect global per-loop metrics
loops <- readRDS("../../output/11JUNE2020_hichip_loopsRaw.mango.rds")
cpm_counts <- cpm(loops@counts)
boo <- (rowMeans(cpm_counts) > 3)
region_annotation <- summary(loops)[boo,c("region")]
sum(boo)

pheatmap::pheatmap(cor(cpm_counts[boo,]))

# Build the per-sample meta data
names <- colnames(cpm_counts)
resp_status = ifelse(grepl("nr", names), "nonresponder", "responder")
celltype = ifelse(grepl("CD4", names), "CD4", "CD8")
timepoint = ifelse(grepl("pre", names), "pre", "post")
phenotype = ifelse(substr(names, nchar(names), nchar(names)) == "M", "memory", "naive")

meta_df <- data.frame(
  resp_status,
  celltype,
  timepoint,
  phenotype
)

library(DESeq2)
library(apeglm)
ddsMat <- DESeqDataSetFromMatrix(countData = loops@counts[boo,],
                                 colData = meta_df,
                                 design = ~ resp_status + celltype + timepoint + phenotype)

ddsMat <- DESeq(ddsMat)

# Supply contrasts to look at how each variable impacts the result
rdf_resp_status <- results(ddsMat, contrast = c("resp_status", "nonresponder", "responder"))
rdf_celltype <- results(ddsMat, contrast = c("celltype", "CD4", "CD8"))
rdf_timepoint <- results(ddsMat, contrast = c("timepoint", "pre", "post"))
rdf_phenotype <- results(ddsMat, contrast = c("phenotype", "memory", "naive"))

sum(rdf_resp_status$padj < 0.1, na.rm = TRUE)
sum(rdf_celltype$padj < 0.1, na.rm = TRUE)
sum(rdf_timepoint$padj < 0.1, na.rm = TRUE)
sum(rdf_phenotype$padj < 0.1, na.rm = TRUE)

dds_shrink <- lfcShrink(ddsMat, coef="resp_status_responder_vs_nonresponder", type="ashr")
plotMA(dds_shrink, ylim = c(-5, 5))


