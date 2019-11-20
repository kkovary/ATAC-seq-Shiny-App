library(shiny)
library(tidyverse)
#library(shinycssloaders)
library(shinysky)
source('Gviz_Plots.R')
library(DT)
library(shinyWidgets)


# Set Bioconductor repositories for shinyapps.io:
# library(BiocManager)
# options(repos = BiocManager::repositories())

# Load data

# load in genelist to speed loading up
# reading csv is faster than rds
data <- readRDS('lite_bc_annotated_rna_dataframe_long.RDS') %>%
  as_tibble()

peaks.gr <- readRDS('All_Merged_Peaks_GenomicRanges.RDS')
cor.gr <- readRDS('correlation_genomic_ranges.RDS')
motifs <- readRDS('Motif_SummarizedExperiment.RDS')

# Functions

plotrna <- function(gene_id, annotated.counts = rna.df.gg, stat.method = "loess", 
                    specify.transcript = NULL, display.se = F, y.limits = c(0,NA), span = 1, log = T) {
  annotated.counts <- annotated.counts %>% filter(gene.symbol==gene_id)
  if (!(is.null(specify.transcript))) {
    annotated.counts <- filter(annotated.counts, transcript_id %in% specify.transcript)
  } 
  if (log==T) {
    gg <- ggplot(annotated.counts, aes(x = Age_day, y = log2.TPM, color = Specific.Type))
  } else {
    gg <- ggplot(annotated.counts, aes(x = Age_day, y = TPM, color = Specific.Type))
  }
  
  gg <- gg +
    geom_point(size = 3) +
    stat_smooth(method = stat.method, aes(fill = Specific.Type), se = display.se, alpha = 0.15, span = span) +
    scale_color_manual(values = c("#8C2D04", "#5f2a99", "#005A32")) +
    scale_fill_manual(values = c("#8C2D04", "#5f2a99", "#005A32")) +
    theme_classic() +
    xlim(c(0,600)) +
    scale_y_continuous(limits = y.limits) +
    ggtitle(label = sprintf("%s RNA expression", unique(annotated.counts$gene.symbol))) +
    theme(aspect.ratio = 1, legend.position="bottom") + 
    ylab('log2(TPM)') + xlab('Age (Day)') +
    facet_wrap( ~ transcript_id, scales = "free_y", ncol = 4)
  return(gg)
}

