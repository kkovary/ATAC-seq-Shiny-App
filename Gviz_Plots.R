###############################################################################
# VISUALIZATION of ATAC data through R
#
#  Gviz exploration
#
#
#
###############################################################################
# Setup
###############################################################################
library(Gviz)
library(GenomicRanges)
library(S4Vectors)
library(SummarizedExperiment)
# colors <- c('#E181F4','#7DCD2C','#F6D7B5','#F9DAFF','#D1E8BA',
#             '#EFE90D','#9404B4','#DC7511','#4C9006','#CB21ED',
#             '#F2A760','#F2B9FF','#B0E57C')


attachMotifs <-
  function(motif.list,
           se,
           chr = NULL,
           beg = NULL,
           END = NULL) {
    se <- se[seqnames(motifs) == chr &
               start(motifs) >= beg &
               end(motifs) <= END]
    
    output <- lapply(motif.list, function(x) {
      granges(se)[assay(se)[, x]]
    })
  }

plotGenomeView <- function(gene.symbol = GENE,
                           slop = SLOP,
                           genome = "hg38",
                           anno.gr = peaks.gr,
                           coverage.list,
                           ylims = c(0, 100),
                           chr = NULL,
                           beg = NULL,
                           END = NULL,
                           transcriptID = NULL,
                           corCut = 0,
                           pval_cut = 1,
                           cluster_id = NULL,
                           motifs_list = NULL,
                           selected_rows = NULL,
                           greater_less = 'Greater than',
                           snp_xrange = NULL) {
  
  print("obtaining coordinates")
  
  if (is.null(beg) & is.null(END) & length(chr) == 0) {
    coords <- getGtfCoords(gene.symbol, ENSEMBL_hg38_local_fromGTF)
    
    chr <- coords[[1]]
    beg <- coords[[2]]
    END <- coords[[3]]
  }
  
  
  axisTrack <- GenomeAxisTrack(fontsize = 20)
  idxTrack <-
    IdeogramTrack(genome = genome,
                  chromosome = chr,
                  fontsize = 20)
  
  print("filtering")
  filt.gr <-
    anno.gr[seqnames(anno.gr) == chr &
              start(anno.gr) > beg - slop &
              end(anno.gr) < END + slop]
  
  print("importing peaks")
  peakTrack <-
    AnnotationTrack(filt.gr,
                    name = "Distal Peaks",
                    col = 'transparent',
                    size = 3)
  
  print("importing coverage")
  covTrackList <- lapply(1:length(coverage.list), function(x) {
    DataTrack(
      range = coverage.list[[x]],
      genome = genome,
      type = "histogram",
      name = color.scheme$name2[x],
      chromosome = chr,
      start = beg,
      end = END,
      col.histogram = colors[x],
      fill.histogram = colors[x],
      col.axis = "black",
      background.title = colors[x],
      fontsize = 12
    )
  })
  
  # Subset the correlation table
  if (length(cor.gr[elementMetadata(cor.gr)$correlated.gene.symbol %in% gene.symbol]) > 0 |
      gene.symbol[1] == 'Any') {
    cor.gr.subset <-
      corFilter(
        cor_table = cor.gr,
        greater_less = greater_less,
        cor_cut = corCut,
        transcript_ID = transcriptID,
        cluster_id = cluster_id,
        pval_cut = pval_cut,
        beg = beg,
        END = END,
        chr = chr,
        gene.name = gene.symbol
      )
    
    corTrack <- AnnotationTrack(
      cor.gr.subset,
      name = 'Cor',
      fill = elementMetadata(cor.gr.subset)$colors,
      col = 'transparent'
    )
  } else{
    corTrack <- NULL
  }
  
  
  if (length(selected_rows) > 0) {
    selected_rows_track <- AnnotationTrack(
      cor.gr.subset[selected_rows,],
      name = 'Selected',
      fill = elementMetadata(cor.gr.subset[selected_rows,])$colors,
      #col = '#DCDCDC'
      col = 'transparent',
      size = 3
    )
  } else{
    selected_rows_track <- NULL
  }
  
  
  
  if (!is.null(motifs_list)) {
    motifs_tracks <- attachMotifs(motifs_list, motifs, chr, beg, END)
    
    motifsTrackList <- lapply(1:length(motifs_tracks), function(x) {
      AnnotationTrack(
        range = motifs_tracks[[x]],
        fill = pal_jco()(10)[x],
        #col = '#DCDCDC',
        col = 'transparent',
        name = motifs_list[x]
      )
    })
  } else{
    motifsTrackList <- NULL
  }
  
  
  plotList <-
    c(
      idxTrack,
      axisTrack,
      covTrackList,
      peakTrack,
      selected_rows_track,
      corTrack,
      ENSEMBL_Gviz_GeneTrackRegionObject,
      motifsTrackList
    )
  
  
  if (gene.symbol[1] == 'Any' & is.null(snp_xrange)) {
    highlight <- data.frame(transcript_start = NULL,
                            transcript_end = NULL)
  } else if (!is.null(snp_xrange)) {
    highlight <- data.frame(transcript_start = snp_xrange$start,
                            transcript_end = snp_xrange$end)
  } else{
    if (length(transcriptID) > 1) {
      highlight <-
        data.frame(
          transcript_start = getGtfCoords(gene.symbol, ENSEMBL_hg38_local_fromGTF)[[2]],
          transcript_end = getGtfCoords(gene.symbol, ENSEMBL_hg38_local_fromGTF)[[3]]
        )
    } else{
      highlight <-
        transcript_locations %>% filter(refseq_mrna == transcriptID)
    }
  }
  
  
  print("plotting")
  plotTracks(
    HighlightTrack(
      plotList,
      start = highlight$transcript_start,
      end = highlight$transcript_end,
      chr = chr,
      col = '#EEEEEE',
      fill = '#EEEEEE'
    ),
    transcriptAnnotation = "gene",
    ylim = ylims,
    col.title = 'black',
    from = beg,
    to = END,
    title.width = 3,
    lwd = 2,
    min.height = 10
  )
  
}


###################################################
# Plot clusters and motifs tracks
###################################################

# plot_clust_motif <- function(gene.symbol = GENE, slop = SLOP,
#                              genome = "hg38", anno.gr = peaks.gr,
#                              coverage.list, ylims = c(0,100),
#                              coords = NULL, chr = NULL,
#                              beg = NULL, END = NULL,
#                              transcriptID = NULL,
#                              corCut = 0,
#                              pval_cut = 1,
#                              cluster_id = NULL,
#                              motifs_list = NULL){
#
#
#   biomTrack <- coords[[4]]
#   filt.gr <- anno.gr[seqnames(anno.gr)==chr & start(anno.gr) > beg - slop & end(anno.gr) < END + slop]
#   peakTrack <- AnnotationTrack(filt.gr, name = "Distal Peaks", col = 'transparent')
#
#   if(is.na(corCut)){
#     cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == transcriptID) &
#                               (elementMetadata(cor.gr)$cluster.name %in% cluster_id) &
#                               (elementMetadata(cor.gr)$vs.null.p.value <= pval_cut) &
#                               start(cor.gr) > beg &
#                               end(cor.gr) < END]
#   } else{
#     cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == transcriptID) &
#                               (elementMetadata(cor.gr)$estimate >= corCut) |
#                               (elementMetadata(cor.gr)$estimate <= -corCut) &
#                               (elementMetadata(cor.gr)$cluster.name %in% cluster_id) &
#                               (elementMetadata(cor.gr)$vs.null.p.value <= pval_cut)&
#                               start(cor.gr) > beg &
#                               end(cor.gr) < END]
#   }
#
#
#   corTrack <- AnnotationTrack(
#     cor.gr.subset,
#     name = 'Cor',
#     fill = elementMetadata(cor.gr.subset)$cluster.color,
#     col = '#DCDCDC'
#     #col = 'transparent'
#   )
#
#   if(!is.null(motifs_list)){
#     motifs_tracks <- attachMotifs(motifs_list, motifs, chr, beg, END)
#
#     motifsTrackList <- lapply(1:length(motifs_tracks), function(x) {
#       AnnotationTrack(range = motifs_tracks[[x]],
#                       fill = pal_jco()(10)[x],
#                       col = '#DCDCDC',
#                       name = motifs_list[x])
#     })
#
#
#
#   }
#
#   if(is.null(motifs_list)){
#     plotList <- c(peakTrack, corTrack, biomTrack)
#   } else{
#     plotList <- c(peakTrack, corTrack, biomTrack, motifsTrackList)
#   }
#
#   highlight <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol",'transcript_start','transcript_end'),
#                      filters = "refseq_mrna", values = transcriptID, mart= ensembl)
#
#   plotTracks(HighlightTrack(plotList,
#                             start = highlight$transcript_start,
#                             end = highlight$transcript_end,
#                             chr = chr),
#              transcriptAnnotation = "name",
#              ylim = ylims, col.title = 'black', from = beg, to = END,
#              title.width = 3)
#
# }

###############################################################################
# MAIN
###############################################################################

# Get BiomaRt objects

color.scheme <- read_csv('Data/Color_Scheme.csv')

coverage.files <- color.scheme$file
strsplits <- str_split(coverage.files, "_")
coverage.names <-
  sapply(strsplits, function(x)
    paste(x[c(5, 2, 3)], collapse = " "))
coverage.list <- as.list(coverage.files)
names(coverage.list) <- coverage.names

colors = color.scheme$color
