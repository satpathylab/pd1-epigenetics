#' ---
#' title: "PD1 chromVAR analysis"
#' author: "Caleb Lareau"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' ### The goal of this document is to produce several interactive plots for ATAC data
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
# Import libraries
library(plotly)
library(heatmaply)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)
library(BuenColors)
library(chromVARxx)
library(chromVARmotifs)
library(data.table)

#' ## Initialize parallel processing
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
library(BiocParallel)
register(MulticoreParam(4)) # adjust according to your machine

#' ## Import/Filter data
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE, eval = TRUE
peakfile <- "../../data/fromMA/peak-ids.bed"
peaks <- getPeaks(peakfile)
counts_ma_in <- fread("../../data/fromMA/combined-counts-clean-subject-corrected-quantiles-ma.tsv")
metadata_ma <- fread("../../data/fromMA/metadata-clean-ma.tsv")
counts_mat <- round(data.matrix(data.frame(counts_ma_in)[,-60]),0)

# Make Summarized Experiment
counts <- SummarizedExperiment(
  rowRanges = peaks, 
  colData = metadata_ma,
  assays = list(counts = counts_mat)
)

#' ## Get GC content/peak; get motifs from chromVARmotifs package
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
counts <- addGCBias(counts, genome = BSgenome.Hsapiens.UCSC.hg19)
data("human_pwms_v2") # also human_pwms_v1 with more motifs (typically redunant)
motif_ix <- matchMotifs(human_pwms_v2, counts, genome = BSgenome.Hsapiens.UCSC.hg19)

#' ## Compute deviations; typically time consuming
#+ cache = TRUE, message=FALSE, warning=FALSE, echo = TRUE
dev <- computeDeviations(object = counts, annotations = motif_ix)

#' ## Find variable motifs
#+ cache = FALSE, message=FALSE, warning=FALSE, echo = TRUE, fig.align='center', fig.cap = "Figure 1a; variability across motifs"
variabilityAll <- computeVariability(dev)
plotVariability(variabilityAll, use_plotly = TRUE)

#' ## Examine variable motifs in sample
#+ cache = FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 20, echo = TRUE, fig.align='center', fig.cap = "Figure 2a; motifs x samples"
mostvariable <- tail(sort(variabilityAll$variability, index.return = TRUE)$ix,30)
m <- assays(dev)[["z"]][mostvariable,]
rownames(m) <- rowData(dev)$name[mostvariable]
heatmaply(m, colors = jdb_palette("solar_extra"))

#' ## Plot in tSNE
#' #+ cache = FALSE, message=FALSE, warning=FALSE, fig.height = 10, fig.width = 10, echo = TRUE, fig.align='center', fig.cap = "Figure 3a; sample tSNE in motif space"
plotDeviationsTsne(dev, deviationsTsne(dev, perplexity = 15, threshold = 2.5), sample_column = "cell.type",
                    shiny = FALSE)[[1]] %>% plotly::ggplotly()

# Export data
saveRDS(round(assays(dev)[["z"]], 2), file = "../../output/27JULY_chromvar_execution.rds")

