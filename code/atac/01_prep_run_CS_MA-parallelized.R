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

# keep_peaks <- !(in_blacklist | near_TSS | in_sex_chromosome)
keep_peaks <- !(in_blacklist | in_sex_chromosome)
table(keep_peaks)

# Import and filter counts to score samples
counts_ma_in <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/combined-counts-clean-ma.tsv")
metadata_ma <- fread("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/metadata-clean-ma.tsv")


counts_ma_in <- counts_ma_in %>%
  dplyr::select(grep("CD4", colnames(.)))
metadata_ma <- metadata_ma %>%
  dplyr::filter(grepl("CD4", sample.id ))
stopifnot(all(metadata_ma[["sample.id"]] == colnames(counts_ma_in)[-60]))

# Drop the peak index and make an actual data matrix object
pd1_samples <- cpm(data.matrix(data.frame(counts_ma_in))[keep_peaks, ])

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
# signature_matrix <- cpm(counts_sp_pop[keep_peaks,c(-3,-4)]) # remove the th2 sample
signature_matrix <- cpm(counts_sp_pop[keep_peaks, c("Naive","Th1","Th2","Th17","Treg")]) # keep Th2 and Th1-17 sample


RUNCIBERSORT <- function(signature_matrix, pd1_samples, variable.peaks) {
  # Subset to variable peaks
  row_vars <- rowVars(signature_matrix)
  cutoff <- sort(row_vars, decreasing = TRUE)[variable.peaks]
  
  signature_matrix_ss <- log1p(signature_matrix[row_vars >= cutoff,])
  pd1_samples_ss <- log1p(pd1_samples[row_vars >= cutoff,])
  
  # Run CIBERSORT on input mixture (pd1_samples) using signature matrix (signature_matrix)
  resmat <- sapply(1:dim(pd1_samples_ss)[2], model_sample_cibersort, pd1_samples_ss, signature_matrix_ss)
  colnames(resmat) <- colnames(pd1_samples_ss)
  resdf <- data.frame(celltype = colnames(signature_matrix_ss), resmat)
  write.table(resdf, paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/18JUNE2020_CIBERSORT_ATAC_", variable.peaks,".tsv", sep= ""),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


###################################################
## setup parallel backend to use many processors ##
###################################################
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

finalMatrix <- foreach(i=c(5000,10000,15000)) %dopar% {
  library(dplyr)
  library(data.table)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(stringr)
  library(matrixStats)
  library(ggplot2)
  library(EnsDb.Hsapiens.v75)
  
  setwd("/Users/mamouzgar/data")
  source("hichip-for-io/meelad/deconvolution/pd1-epigenetics/code/atac/00_CIBERSORT.R")
  
  variable.peaks = i
  tempMatrix = RUNCIBERSORT(signature_matrix, pd1_samples, variable.peaks) #calling a function
  #do other things if you want
  
  # tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
#stop cluster
stopCluster(cl)




# ## begin loop through different cut-offs 
# future.apply::future_lapply(c(5000,10000,15000), function (x) {
#   print(x)
#   variable.peaks <- x
#   # Subset to variable peaks
#   row_vars <- rowVars(signature_matrix)
#   cutoff <- sort(row_vars, decreasing = TRUE)[x]
# 
#   signature_matrix_ss <- log1p(signature_matrix[row_vars >= cutoff,])
#   pd1_samples_ss <- log1p(pd1_samples[row_vars >= cutoff,])
# 
#   # Run CIBERSORT on input mixture (pd1_samples) using signature matrix (signature_matrix)
#   resmat <- sapply(1:dim(pd1_samples_ss)[2], model_sample_cibersort, pd1_samples_ss, signature_matrix_ss)
#   colnames(resmat) <- colnames(pd1_samples_ss)
#   resdf <- data.frame(celltype = colnames(signature_matrix_ss), resmat)
#   write.table(resdf, paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/18JUNE2020_CIBERSORT_ATAC_", variable.peaks,".tsv", sep= ""),
#               sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# })
# 
# for (variable.peaks in c(5000,10000,15000)) { 
# # Subset to variable peaks
# row_vars <- rowVars(signature_matrix)
# cutoff <- sort(row_vars, decreasing = TRUE)[5000]
# 
# signature_matrix_ss <- log1p(signature_matrix[row_vars >= cutoff,])
# pd1_samples_ss <- log1p(pd1_samples[row_vars >= cutoff,])
# 
# # Run CIBERSORT on input mixture (pd1_samples) using signature matrix (signature_matrix)
# resmat <- sapply(1:dim(pd1_samples_ss)[2], model_sample_cibersort, pd1_samples_ss, signature_matrix_ss)
# colnames(resmat) <- colnames(pd1_samples_ss)
# resdf <- data.frame(celltype = colnames(signature_matrix_ss), resmat)
# write.table(resdf, paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/18JUNE2020_CIBERSORT_ATAC_", variable.peaks,".tsv", sep= ""), 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# }
# Visualize results
# melt_df <- reshape2::melt(resdf, id.vars = "celltype")
# df <- data.frame(melt_df, str_split_fixed(melt_df$variable, "_", 3))
# 
# p1 <- ggplot(df, aes(x = as.numeric(as.factor(X1)), y = value, fill = celltype)) +
#   geom_bar(stat = "identity", position = "stack", color = "black") + theme_bw()+
#   facet_grid(X2~X3) + theme(legend.position = "bottom") +
#   # scale_fill_manual(values = c("dodgerblue2", "firebrick", "purple3", "green4")) +
#   labs(x = "Donor #", y = "Estimated celltype proportions") +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(breaks = c(1:8))
# ggsave(p1, file = "hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/17JUNE2020_StackedBarPlot_ATAC.pdf", width = 5, height = 5, useDingbats=FALSE)
