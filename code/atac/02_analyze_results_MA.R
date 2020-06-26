library(tidyverse)
library(magrittr)
library(ggpubr)
setwd("/Users/mamouzgar/data")


output.files <- list.files("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output", full.names = TRUE)
output.files <- output.files[grepl("CIBERSORT", output.files)]
output.files <- output.files[!(grepl("subject|11JUNE",output.files))]


meta <- read.csv("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/metadata-clean-ma.tsv",
                 sep = "\t",
                 check.names = FALSE)

for (filepath in output.files) { 
  
  method = basename(filepath)
  method = gsub(".tsv", "",method)
  print(method)
  resdf <- read.csv(filepath,
                    sep = "\t",
                    check.names = FALSE)
  
  # Visualize results
  melt_df <- resdf %>%
    tidyr::gather(key = "sample.id", value = "percent",-celltype) %>%
    dplyr::filter(grepl("CD4", sample.id)) %>%
    mutate(sample.id = gsub("X", "",sample.id),
           sample.id = gsub("\\.","-",sample.id)) %>%
    left_join(.,meta, by = "sample.id") %>%
    mutate(cell.type.subtype = cell.type,
           cell.type = gsub("N|M", "",cell.type.subtype),
           cell.subtype = case_when(grepl("N",cell.type.subtype) ~ "naive",
                                    grepl("M", cell.type.subtype) ~ "memory"),
           
           ## create factor levels
           cell.type = factor(cell.type, levels = c("CD4")),
           cell.type.subtype = factor(cell.type.subtype, levels = c("CD4N","CD4M")),
           bor = factor(bor, levels = c("POD","PR")),
           timepoint = factor(timepoint, levels = c("pretreatment","posttreatment")),
           
           subject_timepoint = paste(subject, timepoint, sep = "_"))
  
  
  pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", method , ".pdf", sep=""),
      useDingbats = FALSE,
      width = 6,
      height = 6)
  print( ggbarplot(melt_df, x = "subject_timepoint", y = "percent",
                   subtitle = method,
                   color = "celltype",
                   palette = "Dark2",
                   facet.by = "cell.type.subtype") +
           rotate_x_text(90) )
  dev.off()
  
  pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", method , "timepoint-unpaired.pdf", sep=""),
      useDingbats = FALSE,
      width = 6,
      height = 9)
  print( ggpaired(melt_df, x = "timepoint", y = "percent",
                  id = "subject",
                  subtitle = method,
                  line.color = "gray",
                  color = "celltype",
                  palette = "Dark2",
                  facet.by = c("celltype", "cell.type.subtype")) +
           stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.format") )
  dev.off()
  pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", method , "-timepoint-paired.pdf", sep=""),
      useDingbats = FALSE,
      width = 6,
      height = 9)
  print( ggpaired(melt_df, x = "timepoint", y = "percent",
                  id = "subject",
                  subtitle = method,
                  line.color = "gray",
                  color = "celltype",
                  palette = "Dark2",
                  facet.by = c("celltype", "cell.type.subtype")) +
           stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format") )
  dev.off()
  
  pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", method , "CD4N-response-unpaired.pdf", sep=""),
      useDingbats = FALSE,
      width = 6,
      height = 9)
  print ( ggpaired(dplyr::filter(melt_df, cell.type.subtype =="CD4N"), x = "bor", y = "percent",
                   subtitle = method,
                   # id = "subject",
                   line.color = "gray",
                   color = "celltype",
                   palette = "Dark2",
                   facet.by = c("celltype", "timepoint")) +
            stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.format") )
  dev.off()
  
  pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", method , "CD4M-response-unpaired.pdf", sep=""),
      useDingbats = FALSE,
      width = 6,
      height = 9)
  print ( ggpaired(dplyr::filter(melt_df, cell.type.subtype =="CD4M"), x = "bor", y = "percent",
                   subtitle = method,
                   # id = "subject",
                   line.color = "gray",
                   color = "celltype",
                   palette = "Dark2",
                   facet.by = c("celltype", "timepoint")) +
            stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.format") )
  dev.off()

  
  
  
}



#########
## PCA ##
#########








# resdf <- read.csv("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/17JUNE2020_CIBERSORT_ATAC.tsv",
#             sep = "\t",
#             check.names = FALSE)
# 
# meta <- read.csv("hichip-for-io/meelad/deconvolution/pd1-epigenetics/data/fromMA/metadata-clean-ma.tsv",
#                   sep = "\t",
#                   check.names = FALSE)
# 
# # Visualize results
# melt_df <- resdf %>%
#   tidyr::gather(key = "sample.id", value = "percent",-celltype) %>%
#   mutate(sample.id = gsub("X", "",sample.id),
#          sample.id = gsub("\\.","-",sample.id)) %>%
#   left_join(.,meta, by = "sample.id") %>%
#   mutate(cell.type.subtype = cell.type,
#          cell.type = gsub("N|M", "",cell.type.subtype),
#          cell.subtype = case_when(grepl("N",cell.type.subtype) ~ "naive",
#                                   grepl("M", cell.type.subtype) ~ "memory"),
# 
#          ## create factor levels
#          cell.type = factor(cell.type, levels = c("CD4")),
#          cell.type.subtype = factor(cell.type.subtype, levels = c("CD4N","CD4M")),
#          bor = factor(bor, levels = c("POD","PR")),
#          timepoint = factor(timepoint, levels = c("pretreatment","posttreatment")),
# 
#          subject_timepoint = paste(subject, timepoint, sep = "_"))
# 
# 
# pdf(paste("hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/deconvolution-plots/", Sys.Date(), ,sep=""),
#     useDingbats = FALSE,
#     width = 6,
#     height = 6)
# ggbarplot(melt_df, x = "subject_timepoint", y = "percent",
#           subtitle =
#           color = "celltype",
#           palette = "Dark2",
#           facet.by = "cell.type.subtype") +
#   rotate_x_text(90)
# dev.off()
# 
# ggpaired(melt_df, x = "timepoint", y = "percent",
#          id = "subject",
#          line.color = "gray",
#          color = "celltype",
#          palette = "Dark2",
#          facet.by = c("celltype", "cell.type.subtype")) +
#   stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format")
# 
# ggpaired(dplyr::filter(melt_df, cell.type.subtype =="CD4M"), x = "bor", y = "percent",
#          id = "subject",
#          line.color = "gray",
#          color = "celltype",
#          palette = "Dark2",
#          facet.by = c("celltype", "timepoint")) +
#   stat_compare_means(method = "wilcox.test", paired = FALSE, label = "p.format")
# 
# 
# ggsave(p1, file = "hichip-for-io/meelad/deconvolution/pd1-epigenetics/output/17JUNE2020_StackedBarPlot_ATAC.pdf", width = 5, height = 5, useDingbats=FALSE)
