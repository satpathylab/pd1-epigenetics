library(data.table)
library(Matrix)
library(GenomicRanges)
library(Hmisc)
library(data.table)
library(zoo)
library(dplyr)
library(plotrix)
library(BuenColors)
library(stringr)

options(scipen=999)

# Specify the directory

# CL Mac
base_dir <- "/Users/clareau/sherlock_oak/users/mamouzgar/HiChIP_Pipeline_local/io-chromatin/output/"

# Misc. functions
importTable <- function(file){
  print(file)
  return(data.frame(fread(file, sep = "\t", header =  FALSE)))
}

# Pick your smoothing window (1 means no smoothing, 2 averages over 2 adjacent bins, etc. )
rollmean_win <- 3
resolution <- 5000 # for what we get from HiC-Pro

# Import coordinate for one particular sample as this should be all the same given some resolution
bed_df <- importTable(paste0(base_dir,"/hichip.rpost.CD4M/matrices/hichip.rpost.CD4M_",as.character(resolution),"_abs.bed"))
colnames(bed_df) <- c("chr", "start", "end", "idx")
bed_gr <- makeGRangesFromDataFrame(bed_df, keep.extra.columns = TRUE) 
coord_vec <- bed_df$start; names(coord_vec) <- as.character(bed_df$idx)

# Import matrix tables and then put them as a list of arbitrary length
samples <- c("nrpost.CD4M","nrpost.CD4N","nrpost.CD8M","nrpost.CD8N",
             "nrpre.CD4M","nrpre.CD4N","nrpre.CD8M","nrpre.CD8N",
             "rpost.CD4M","rpost.CD4N","rpost.CD8M","rpost.CD8N",
             "rpre.CD4M","rpre.CD4N","rpre.CD8M","rpre.CD8N")

matrix_files <- paste0(base_dir, "/hichip.", samples, "/matrices/hichip.", samples,"_", resolution, ".matrix")

# Import sparse matrix counts 
samples_list <- lapply(matrix_files, importTable); names(samples_list) <- samples
saveRDS(samples_list, file = "../../../hichip_meelad_1kb_res.rds")

# Function takes a GRanges of the v4C view point and the width around the viewpoint (to be padded on both sides)
make_v4C_plot <- function(v4C_one, width_around_viewpoint= 200000){
  
  # Specify one region from candidates 
  ov <- findOverlaps(v4C_one, bed_gr, minoverlap = width(v4C_one)/2)
  hicpro_idx  <- subjectHits(ov)
  chromosome <- as.character(seqnames(v4C_one))
  
  # Function to grab the interaction counts per replicate 
  processTable <- function(table, center_coord, sampleID){
    
    # Do a normalization
    table_norm <- sum(table$V3)/1000000 # effectively a CPM normalization
    counts1_left=table[table[,2] == center_coord,c(1,3)]
    counts1_right=table[table[,1] == center_coord,c(2,3)]
    colnames(counts1_left) = colnames(counts1_right) = c("idx","signal")
    
    # Combine these and filter for a chromosome
    counts1_merge=rbind(counts1_left, counts1_right)
    counts1_merge_chr <- counts1_merge[counts1_merge[,1] %in% bed_df[bed_df$chr == chromosome,"idx"],]
    counts1_merge_chr$bp <- coord_vec[as.character(counts1_merge_chr$idx)]
    counts1_merge <- counts1_merge_chr[,c("bp", "signal")]
    counts1_merge_ord = counts1_merge[order(counts1_merge[,1]),]
    counts1_merge_ord[,2] = counts1_merge_ord[,2]/table_norm
    counts1_merge_ord[,1] = counts1_merge_ord[,1] + resolution/2
    counts1_merge_ord$line = suppressWarnings(rollmean(counts1_merge_ord[,2], rollmean_win, fill=c(0,0,0)))
    counts1_merge_ord$SampleID <-sampleID
    return(counts1_merge_ord)
  }
  
  # Apply over all samples to create a long interaction data frame
  interactDf <- lapply(names(samples_list), function(sample_name){
    processTable(samples_list[[sample_name]], hicpro_idx, sample_name)
  }) %>% rbindlist() %>% data.frame()
  
  # Specify the plot conditions
  xlim_minus=start(v4C_one) - width_around_viewpoint # bp upstream of anchor
  xlim_plus=end(v4C_one) + width_around_viewpoint  # bp downstream of anchor
  gene <- as.character(mcols(v4C_one)$gene)
  
  # Set up plotting stuff
  midpoint <- mean(c(start(v4C_one),end(v4C_one)))
  xlim_par= c(xlim_minus,xlim_plus)
  ymax <- max(interactDf$line)
  ylim_par=c(-0.001,ymax) # y axis limits
  lwd_par = 2 # line thickness
  xlab=paste0("Position on ", chromosome)
  ylab="v4C normalized score"
  plot_name=paste0(gene, " v4C; width: ", round((xlim_plus-xlim_minus)/1000), "kb")
  
  # A few additional meta data variables
  two <- str_split_fixed(interactDf$SampleID, "[.]",2)
  interactDf$celltype <- two[,2]
  interactDf$biotype <- two[,1]
  
  # ggplot implementation
  p1 <- ggplot(interactDf, aes(x = bp, y = line, color = biotype)) + 
    geom_line() + pretty_plot() + 
    scale_color_manual(values = c("firebrick", "orange2", "dodgerblue4", "dodgerblue")) +
    geom_vline(xintercept = midpoint, linetype="dashed") +
    theme(legend.position = "bottom") +
    facet_wrap(~celltype, ncol = 2) +
    coord_cartesian(ylim = ylim_par, xlim = xlim_par) +
    labs(x = xlab, y = ylab, color = "") + ggtitle(plot_name)

  cowplot::ggsave2(p1, file= paste0("../../output/v4c_plots/", gene, "-v4C_hichip.pdf"), height = 7, width = 14)
  gene
}

# Import candidates - bed file of the following:
# chr, start, end, gene/whatever name for the view points
v4c <-data.frame(
  chr = c("chr1", "chr16", "chr12", "chr2", "chr2", "chr6"), 
  start = c(59248739, 50402686, 6909299, 86990371, 192016199, 91006580),
  end = c(59248740, 50402687, 6909300, 86990372, 192016200, 91006581),
  gene = c("JUN", "BRD7", "CD4", "CD8A", "STAT4", "BACH2")
)

v4c_gr <- makeGRangesFromDataFrame(v4c, keep.extra.columns = TRUE)
lapply(1:dim(v4c)[1], function(i){
  make_v4C_plot(v4c_gr[i])
})