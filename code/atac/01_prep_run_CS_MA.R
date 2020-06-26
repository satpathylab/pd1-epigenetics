library(dplyr)
library(data.table)
library(rtracklayer)
library(GenomeInfoDb)
library(stringr)
library(matrixStats)
library(ggplot2)
library(EnsDb.Hsapiens.v75) # hg19 based
setwd("/Users/mamouzgar/data")
source("hichip-for-io/meelad/deconvolution/pd1-epigenetics/code/atac/00_CIBERSORT.R")

cpm <- function(mat){
  t(t(mat)/colSums(mat))*1000000
}

# First, filter the input matrix
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5042844/ under Cibersort methods

# Import peaks + blacklisted peaks + TSS
peaks <- import("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/peak-ids.bed")
bl <- import("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/hg19.full.blacklist.bed")
tss <- promoters(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding", upstream = 1000, downstream = 1000)
seqlevelsStyle(tss) <- 'UCSC'

# Find overlapping peaks
in_blacklist <- 1:length(peaks) %in% queryHits(findOverlaps(peaks,bl))
near_TSS <- 1:length(peaks) %in% queryHits(findOverlaps(peaks,tss))
in_sex_chromosome <- as.logical(seqnames(peaks) %in% c("chrX", "chrY"))

keep_peaks <- !(in_blacklist | near_TSS | in_sex_chromosome)
table(keep_peaks)

# Import and filter counts to score samples
counts_ma_in <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/combined-counts-clean-ma.tsv")
metadata_ma <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/metadata-clean-ma.tsv")


counts_ma_in <- dplyr::select(counts_ma_in, grep("CD4", colnames(counts_ma_in)))
metadata_ma <- dplyr::filter(metadata_ma, grepl("CD4", sample.id))
all.equal(metadata_ma$sample.id, colnames(counts_ma_in))

# Perform batch-correction on depth normalized data
depth_norm_data <- cpm(counts_ma_in)
design = model.matrix(~0 + timepoint  +cell.type, metadata_ma)
counts_ma_in_noBatach <- removeBatchEffect(depth_norm_data,
                                           design= design,
                                           batch=metadata_ma$subject)

stopifnot(all(metadata_ma[["sample.id"]] == colnames(counts_ma_in)))

# Make an actual data matrix object
pd1_samples <- data.matrix(counts_ma_in_noBatach[keep_peaks, ]) 

# Import and process the sorted subsets
counts_sp_in <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/Satpathy2019NatMed/2020June10_SatpathyNatMed_counts.tsv") %>% data.frame() %>% data.matrix()
metadata_sp <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/Satpathy2019NatMed/2020June10_SatpathyNatMed_metadata.tsv", header = FALSE)
metadata_sp_all <- data.frame(
  SRR = metadata_sp[[1]],
  str_split_fixed(metadata_sp[[2]], "_", 3)
)

# Sum over replicates
counts_sp_pop <- sapply(unique(metadata_sp_all$X1), function(ss){
  rowSums(counts_sp_in[,metadata_sp_all$X1 == ss])
})
colnames(counts_sp_pop) <- unique(metadata_sp_all$X1)
signature_matrix <- cpm(counts_sp_pop[keep_peaks, ]) # keep Th2 and Th1-17 sample

# Perform feature scaling on both input mixture (pd1_samples) using signature matrix (signature_matrix)
pd1_samples <- scale(pd1_samples, center = FALSE, scale = TRUE)
signature_matrix <- scale(signature_matrix, center = FALSE, scale = TRUE)

# Subset to variable peaks
row_vars <- rowVars(signature_matrix)
cutoff <- sort(row_vars, decreasing = TRUE)[5000]

signature_matrix_ss <- signature_matrix[row_vars >= cutoff,]
pd1_samples_ss <- pd1_samples[row_vars >= cutoff,]

# Run CIBERSORT on input mixture (pd1_samples) using signature matrix (signature_matrix)
resmat <- sapply(1:dim(pd1_samples_ss)[2], model_sample_cibersort, pd1_samples_ss, signature_matrix_ss)
colnames(resmat) <- colnames(pd1_samples_ss)
resdf <- data.frame(celltype = colnames(signature_matrix_ss), resmat)
write.table(resdf, "hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/17JUNE2020_CIBERSORT_ATAC_5000-subject-corrected.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# #########
# #########
# #########
# model_sample_cibersort <- function(idx, Y, sigMat){
#   print(idx)
#   
#   # Get best k
#   mses <- getMSE(Y[,idx,drop = FALSE], sigMat)
#   best_k <- as.numeric(names(which(min(mses) == mses)))
#   
#   # Run model for best k; get coefficients
#   tunedModel <- svm(sigMat, Y[,idx,drop = FALSE], scale = FALSE, kernel = "linear", type = "nu-regression", nu = best_k)
#   
#   # Estimate coefficients; remove negative values; normalize contributions
#   coefTuned <- t(tunedModel$coefs) %*% tunedModel$SV
#   coefTuned[coefTuned < 0 ] <- 0
#   coefTunedNormed <- round(coefTuned /sum(coefTuned), 2)
#   return(tunedModel)
# }
# pd1_samples_ss[ , 1,drop=FALSE]
# test <- model_sample_cibersort(1, pd1_samples_ss, signature_matrix_ss)
# lol <- t(test$coefs) %*% test$SV
# lol[lol<0] <- 0
# #########
# #########
# #########




# Visualize results
melt_df <- reshape2::melt(resdf, id.vars = "celltype")
df <- data.frame(melt_df, str_split_fixed(melt_df$variable, "_", 3))

p1 <- ggplot(df, aes(x = as.numeric(as.factor(X1)), y = value, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "black") + theme_bw()+
  facet_grid(X2~X3) + theme(legend.position = "bottom") +
  # scale_fill_manual(values = c("dodgerblue2", "firebrick", "purple3", "green4")) +
  labs(x = "Donor #", y = "Estimated celltype proportions") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(1:8))
ggsave(p1, file = "hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/17JUNE2020_StackedBarPlot_ATAC_15000.pdf", width = 5, height = 5, useDingbats=FALSE)
