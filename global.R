library(shiny)
library(tidyverse)
library(shinycssloaders)
library(shinysky)
library(DT)
library(shinyWidgets)
library(ggsci)
library(ggpubr)
library(feather)
source('Gviz_Plots.R')

my_names <- c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS')
my_selected <- c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS')
my_colors <- c('#9E0142','#D53E4F','#F46D43','#FFFFBF','#E6F598','#6BAED6',
               '#3288BD','#5E4FA2','#ABDDA4','#66C2A5','#969696')

SLOP = 50000
# Set Bioconductor repositories for shinyapps.io:
# library(BiocManager)
# options(repos = BiocManager::repositories())

# Load data

# load in genelist to speed loading up
# reading csv is faster than rds
data <- read_feather('lite_bc_annotated_rna_dataframe_long.feather')

peaks.gr <- readRDS('All_Merged_Peaks_GenomicRanges.RDS')

cor.gr <- read_feather('correlation_genomic_ranges.feather') %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

motifs <- readRDS('Motif_SummarizedExperiment.RDS')

ENSEMBL_hg38_local_fromGTF <- read_feather('ENSEMBL_hg38_local_fromGTF.feather') %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

ENSEMBL_Gviz_GeneTrackRegionObject <- readRDS('ENSEMBL_Gviz_GeneTrackRegionObject.RDS')

transcript_locations <- read_feather('transcript_locations.feather')

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

my_checkboxGroupInput <- function(variable, label, choices, selected, colors){
  choices_names <- choices
  if(length(names(choices))>0) my_names <- names(choices)
  div(id=variable,class="form-group shiny-input-checkboxgroup shiny-input-container shiny-bound-input",
      HTML(paste0('<label class="control-label" for="',variable,'">',label,'</label>')),
      div( class="shiny-options-group",
           HTML(paste0('<div class="checkbox" style="background-color:', colors,'">',
                       '<label>',
                       '<input type="checkbox" name="', variable, 
                       '" value="', choices, 
                       '"', ifelse(choices %in% selected, 'checked="checked"', ''), 
                       '/>',
                       '<span>', choices_names,'</span>',
                       '</label>',
                       '</div>', collapse = " "))
      )
  )
}