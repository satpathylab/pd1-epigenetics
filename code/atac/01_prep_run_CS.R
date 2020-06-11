library(dplyr)
library(data.table)
library(rtracklayer)
library(GenomeInfoDb)
library(stringr)
library(matrixStats)
library(ggplot2)
library(EnsDb.Hsapiens.v75) # hg19 based

source("00_CIBERSORT.R")

cpm <- function(mat){
  t(t(mat)/colSums(mat))*1000000
}

# First, filter the input matrix
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5042844/ under Cibersort methods

# Import peaks + blacklisted peaks + TSS
peaks <- import("../../data/fromMA/peak-ids.bed")
bl <- import("../../data/hg19.full.blacklist.bed")
tss <- promoters(EnsDb.Hsapiens.v75, filter = ~ gene_biotype == "protein_coding", upstream = 1000, downstream = 1000)
seqlevelsStyle(tss) <- 'UCSC'

# Find overlapping peaks
in_blacklist <- 1:length(peaks) %in% queryHits(findOverlaps(peaks,bl))
near_TSS <- 1:length(peaks) %in% queryHits(findOverlaps(peaks,tss))
in_sex_chromosome <- as.logical(seqnames(peaks) %in% c("chrX", "chrY"))

keep_peaks <- !(in_blacklist | near_TSS | in_sex_chromosome)
table(keep_peaks)

# Import and filter counts to score samples
counts_ma_in <- fread("../../data/fromMA/combined-counts-clean-ma.tsv")
metadata_ma <- fread("../../data/fromMA/metadata-clean-ma.tsv")

stopifnot(all(metadata_ma[["sample.id"]] == colnames(counts_ma_in)[-60]))

# Drop the peak index and make an actual data matrix object
pd1_samples <- cpm(data.matrix(data.frame(counts_ma_in))[keep_peaks,-60])

# Import and process the sorted subsets
counts_sp_in <- fread("../../data/Satpathy2019NatMed/2020June10_SatpathyNatMed_counts.tsv") %>% data.frame() %>% data.matrix()
metadata_sp <- fread("../../data/Satpathy2019NatMed/2020June10_SatpathyNatMed_metadata.tsv", header = FALSE)
metadata_sp_all <- data.frame(
  SRR = metadata_sp[[1]],
  str_split_fixed(metadata_sp[[2]], "_", 3)
)

# Sum over replicates
counts_sp_pop <- sapply(unique(metadata_sp_all$X1), function(ss){
  rowSums(counts_sp_in[,metadata_sp_all$X1 == ss])
})
colnames(counts_sp_pop) <- unique(metadata_sp_all$X1)
signature_matrix <- cpm(counts_sp_pop[keep_peaks,c(-3,-4)]) # remove the th2 sample

# Subset to variable peaks
row_vars <- rowVars(signature_matrix)
cutoff <- sort(row_vars, decreasing = TRUE)[5000]

signature_matrix_ss <- log1p(signature_matrix[row_vars >= cutoff,])
pd1_samples_ss <- log1p(pd1_samples[row_vars >= cutoff,])

# Run CIBERSORT on input mixture (pd1_samples) using signature matrix (signature_matrix)
resmat <- sapply(1:dim(pd1_samples_ss)[2], model_sample_cibersort, pd1_samples_ss, signature_matrix_ss)
colnames(resmat) <- colnames(pd1_samples_ss)
resdf <- data.frame(celltype = colnames(signature_matrix_ss), resmat)
write.table(resdf, "../../output/11JUNE2020_CIBERSORT_ATAC.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Visualize results
melt_df <- reshape2::melt(resdf, id.vars = "celltype")
df <- data.frame(melt_df, str_split_fixed(melt_df$variable, "_", 3))

p1 <- ggplot(df, aes(x = as.numeric(as.factor(X1)), y = value, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "black") + theme_bw()+
  facet_grid(X2~X3) + theme(legend.position = "bottom") +
  scale_fill_manual(values = c("dodgerblue2", "firebrick", "purple3", "green4")) +
  labs(x = "Donor #", y = "Estimated celltype proportions") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(1:8))
ggsave(p1, file = "../../output/11JUNE2020_StackedBarPlot_ATAC.pdf", width = 5, height = 5, useDingbats=FALSE)
