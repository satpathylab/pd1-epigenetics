library(diffloop)
library(matrixStats)
library(dplyr)

cpm <- function(mat){
  t(t(mat)/colSums(mat))*1000000
}

# Collect global per-loop metrics
loops <- readRDS("../../output/11JUNE2020_hichip_loopsRaw.mango.rds")
cpm_counts <- cpm(loops@counts)
loops@rowData$mean <- rowMeans(cpm_counts)
loops@rowData$variance <- rowVars(cpm_counts)
loops@rowData$VMR <- loops@rowData$variance/loops@rowData$mean 

# Create summary dataframe
colnames(cpm_counts) <- paste0(colnames(cpm_counts), "_cpm")
loops@rowData$nr_minus_resp <- rowMeans(cpm_counts[,grepl("nr", colnames(cpm_counts))]) - rowMeans(cpm_counts[,!grepl("nr", colnames(cpm_counts))])
loops_df <- cbind(summary(loops), cpm_counts)

# More targeted differences
loops_df$r_nrCD4N <- loops_df$hichip.rpre.CD4N_cpm - loops_df$hichip.nrpre.CD4N_cpm

loops_df %>% arrange(desc(variance)) %>% head()
