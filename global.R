# Set Bioconductor repositories for shinyapps.io:
# library(BiocManager)
# options(repos = BiocManager::repositories())

# Load data

# load in genelist to speed loading up
# reading csv is faster than rds

library(shiny)
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(shinycssloaders)
library(shinyjs)
library(shinysky)
library(DT)
library(shinyWidgets)
library(feather)
library(ggsci)
library(ggpubr)
source('Gviz_Plots.R')

# Load Global Data
#gene_names <- readRDS('Data/gene_names.RDS')
gene_names <- as.vector(read_feather('Data/gene_names.feather')$gene_names)
chromosomes <- c(paste0('chr', c(1:22, 'X', 'Y')))
motifs <- readRDS('Data/Motif_SummarizedExperiment.RDS')

ENSEMBL_hg38_local_fromGTF <-
  read_feather('Data/ENSEMBL_hg38_local_fromGTF.feather') %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

transcript_locations <-
  read_feather('Data/transcript_locations.feather')

ideo_data <- read_feather('Data/ideo_data.feather')



my_names <- c(
  'Mature Glia 1',
  'Mature Glia 2',
  'Glial Progenitors',
  'Pluripotency',
  'Early Differentiation 1',
  'Early Differentiation 2',
  'Pallial Neuron 1',
  'Pallial Neuron 2',
  'Late Neuron 1',
  'Late Neuron 2',
  'Constitutive'
)
my_selected <-
  c(
    'Mature Glia 1',
    'Mature Glia 2',
    'Glial Progenitors',
    'Pluripotency',
    'Early Differentiation 1',
    'Early Differentiation 2',
    'Pallial Neuron 1',
    'Pallial Neuron 2',
    'Late Neuron 1',
    'Late Neuron 2',
    'Constitutive'
  )
my_colors <-
  c(
    '#9E0142',
    '#D53E4F',
    '#F46D43',
    '#FFFFBF',
    '#E6F598',
    '#6BAED6',
    '#3288BD',
    '#5E4FA2',
    '#ABDDA4',
    '#66C2A5',
    '#969696'
  )

SLOP = 50000

# Functions

plotrna <-
  function(gene_id,
           annotated.counts = rna.df.gg,
           stat.method = "loess",
           display.se = F,
           y.limits = c(0, NA),
           span = 1) {
    annotated.counts <-
      annotated.counts %>% filter(gene.symbol == gene_id)
    
    if (nrow(annotated.counts) > 0) {
      gg <-
        ggplot(annotated.counts,
               aes(x = Age_day, y = log2.TPM, color = Specific.Type)) +
        geom_point(size = 3) +
        stat_smooth(
          method = stat.method,
          aes(fill = Specific.Type),
          se = display.se,
          alpha = 0.15,
          span = span
        ) +
        scale_color_manual(values = c("#8C2D04", "#5f2a99", "#005A32")) +
        scale_fill_manual(values = c("#8C2D04", "#5f2a99", "#005A32")) +
        theme_classic() +
        xlim(c(0, 600)) +
        scale_y_continuous(limits = y.limits) +
        ggtitle(label = sprintf("%s RNA expression", unique(annotated.counts$gene.symbol))) +
        theme(aspect.ratio = 1, legend.position = "bottom") +
        ylab('log2(TPM)') + xlab('Age (Day)') +
        facet_wrap( ~ transcript_id, scales = "free_y", ncol = 4)
    } else{
      gg <-
        ggplot(data.frame()) + ggtitle('No RNA-seq Data Available') + theme_bw()
    }
    
    return(gg)
  }

my_checkboxGroupInput <-
  function(variable, label, choices, selected, colors) {
    choices_names <- choices
    if (length(names(choices)) > 0)
      my_names <- names(choices)
    div(
      id = variable,
      class = "form-group shiny-input-checkboxgroup shiny-input-container shiny-bound-input",
      HTML(
        paste0(
          '<label class="control-label" for="',
          variable,
          '">',
          label,
          '</label>'
        )
      ),
      div(class = "shiny-options-group",
          HTML(
            paste0(
              '<div class="checkbox" style="background-color:',
              colors,
              '">',
              '<label>',
              '<input type="checkbox" name="',
              variable,
              '" value="',
              choices,
              '"',
              ifelse(choices %in% selected, 'checked="checked"', ''),
              '/>',
              '<span>',
              choices_names,
              '</span>',
              '</label>',
              '</div>',
              collapse = " "
            )
          ))
    )
  }


# Input SNP ID and get back list of genes within range (example rs144861725)
snp_gene <-
  function(snp_id,
           table = ENSEMBL_hg38_local_fromGTF,
           snp_tbl = snp_table) {
    positions <- dplyr::filter(snp_tbl, SNP.ID == snp_id)
    start <- positions$start - SLOP
    end <- positions$end + SLOP
    
    table <-
      table[seqnames(table) == as.character(positions$chrom) &
              start(table) > positions$start - SLOP &
              end(table) < positions$end + SLOP &
              elementMetadata(table)$gene_biotype %in% c('protein_coding', 'lncRNA')]
    table <- c('Any', unique(table$gene_name))
    
    return(table)
  }

# Find coordinates of gene
getGtfCoords <- function(GENE, tbl = ENSEMBL_hg38_local_fromGTF) {
  gtf2 <- tbl[tbl@elementMetadata$gene_name == GENE]
  
  g.chr <- min(as.character(seqnames(gtf2)))
  g.start <- min(start(gtf2))
  g.end <- max(end(gtf2))
  
  return(list(g.chr, g.start, g.end))
}

# Find genes in xrange
genes_xrange <-
  function(xmin = NULL,
           xmax = NULL,
           chr = NULL,
           table = ENSEMBL_hg38_local_fromGTF) {
    table <-
      table[seqnames(table) == chr &
              start(table) >= xmin & end(table) <= xmax &
              elementMetadata(table)$gene_biotype %in% c('protein_coding', 'lncRNA')]
    return(c('Any', unique(table$gene_name)))
  }


# Filter cor table
corFilter <-
  function(cor_table = cor.gr,
           greater_less = 'Greater than',
           cor_cut = NA,
           transcript_ID = 'Any',
           cluster_id = my_names,
           pval_cut = 1,
           beg = NULL,
           END = NULL,
           chr = NULL,
           gene.name = 'Any') {
    # Filter based on transcript or show all if "Any" gene is displayed
    if (gene.name[1] != 'Any') {
      cor_table = cor_table[elementMetadata(cor_table)$correlated.transcript.ID %in% transcript_ID]
    } else{
      cor_table = cor_table[elementMetadata(cor_table)$correlated.gene.symbol %in% gene.name]
    }
    
    
    # Filter based on correlation
    if (!is.na(cor_cut)) {
      cor_table <-
        cor_table[!is.na(elementMetadata(cor_table)$Pearson.r.peak.transcript)]
      cor_table = switch(greater_less,
                         'Greater than' = cor_table[elementMetadata(cor_table)$Pearson.r.peak.transcript > cor_cut],
                         'Less than' = cor_table[elementMetadata(cor_table)$Pearson.r.peak.transcript < cor_cut])
    }
    
    # Filter based on cluster ID
    if (length(cluster_id) > 0) {
      cor_table = cor_table[elementMetadata(cor_table)$KM.cluster.name %in% cluster_id]
    }
    
    # Filter based on pvalue
    if (!is.na(pval_cut)) {
      cor_table = cor_table[elementMetadata(cor_table)$correlation.global.model.P <= pval_cut]
    }
    
    # Filter based on beginning
    if (!is.na(beg)) {
      cor_table = cor_table[start(cor_table) > beg]
    }
    
    if (!is.na(END)) {
      cor_table = cor_table[end(cor_table) < END]
    }
    
    return(cor_table)
    
  }
