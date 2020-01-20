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
#library(GenomicFeatures)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(biomaRt)
library(S4Vectors)
library(SummarizedExperiment)
# colors <- c('#E181F4','#7DCD2C','#F6D7B5','#F9DAFF','#D1E8BA',
#             '#EFE90D','#9404B4','#DC7511','#4C9006','#CB21ED',
#             '#F2A760','#F2B9FF','#B0E57C')

###############################################################################
# Input
###############################################################################

# GENE = "SYN1"
# SLOP = 50000
# GENOME = "hg38"

###############################################################################
# Functions
###############################################################################

# getEnsemblCoords <- function(GENE) {
#   # Too slow, don't use
#   bm <- getBM(attributes = c("transcript_start", "transcript_end"), 
#               filters = c("hgnc_symbol"), values = c(GENE), mart = ensembl)
#   g.start <- min(bm$transcript_start)
#   g.end <- max(bm$transcript_end)
#   return(c(g.start,g.end))
# }

getGtfCoords <- function(GENE, tbl = gtf) {
  
  gtf2 <- tbl[tbl@elementMetadata$gene_name==GENE]
  
  g.chr <- min(as.character(seqnames(gtf2)))
  g.start <- min(start(gtf2))
  g.end <- max(end(gtf2))
  
  return(list(g.chr, g.start, g.end))
}

# getGvizCoords <- function(gene.symbol) {
#   
#   biomTrack <- BiomartGeneRegionTrack(genome="hg38", 
#                                       name="ENSEMBL", 
#                                       biomart=ensembl, 
#                                       symbol = gene.symbol)
#   
#   g.chr <- Gviz::seqlevels(biomTrack)
#   pos <- Gviz::position(biomTrack)
#   g.beg <- min(pos)
#   g.end <- max(pos)
#   
#   return(list(g.chr, g.beg, g.end, biomTrack))
# }


attachMotifs <- function(motif.list, se, chr = NULL, beg = NULL, END = NULL) {
  se <- se[seqnames(motifs) == chr & 
             start(motifs) >= beg & 
             end(motifs) <= END]
  
  output <- lapply(motif.list, function(x) {
    # temp <- granges(se)[assay(se)[, x]]
    # temp[start(temp) > beg & end(temp) < END]
    granges(se)[assay(se)[, x]]
  })
}

plotGenomeView <- function(gene.symbol = GENE, slop = SLOP, 
                           genome = "hg38", anno.gr = peaks.gr, 
                           coverage.list, ylims = c(0,100),
                           coords = NULL, chr = NULL,
                           beg = NULL, END = NULL,
                           transcriptID = NULL,
                           corCut = 0,
                           pval_cut = 1,
                           cluster_id = NULL,
                           motifs_list = NULL,
                           selected_rows = NULL) {
  
  #print("creating track")
  #biomTrack <- BiomartGeneRegionTrack(genome=genome, name="ENSEMBL", biomart=ensembl, symbol = gene.symbol)
  
  print("obtaining coordinates")
  if(is.null(coords)){
    coords <- getGtfCoords(gene.symbol, ENSEMBL_hg38_local_fromGTF)
  }
  
  if(is.null(beg) & is.null(END)){
    chr <- coords[[1]]
    beg <- coords[[2]]
    END <- coords[[3]]
  }
  
  
  #biomTrack <- BiomartGeneRegionTrack(genome=genome, name="ENSEMBL", biomart=ensembl, symbol = gene.symbol)
  print(as.character(coords[1:3]))
  
  axisTrack <- GenomeAxisTrack(fontsize = 20)
  idxTrack <- IdeogramTrack(genome = genome, chromosome = chr)
  
  print("filtering")
  filt.gr <- anno.gr[seqnames(anno.gr)==chr & start(anno.gr) > beg - slop & end(anno.gr) < END + slop] 
  
  print("importing peaks")
  peakTrack <- AnnotationTrack(filt.gr, name = "Distal Peaks", col = 'transparent')
  
  print("importing coverage")
  covTrackList <- lapply(1:length(coverage.list), function(x) {
    DataTrack(range = coverage.list[[x]], genome = genome, 
              type = "histogram", name = color.scheme$name2[x], 
              chromosome = chr, start = beg, end = END,
              col.histogram = colors[x],
              fill.histogram = colors[x],
              col.axis = "black",
              background.title = colors[x],
              fontsize = 12)
  })
  
  if(is.na(corCut)){
    cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$correlated.transcript.ID %in% transcriptID) &
                              (elementMetadata(cor.gr)$KM.cluster.name %in% cluster_id) &
                              (elementMetadata(cor.gr)$correlation.local.model.P <= pval_cut) &
                              start(cor.gr) > beg & 
                              end(cor.gr) < END]
  } else{
    cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$correlated.transcript.ID %in% transcriptID) &
                              (elementMetadata(cor.gr)$Pearson.r.peak.transcript >= corCut) |
                              (elementMetadata(cor.gr)$Pearson.r.peak.transcript <= -corCut) & 
                              (elementMetadata(cor.gr)$KM.cluster.name %in% cluster_id) &
                              (elementMetadata(cor.gr)$correlation.local.model.P <= pval_cut)&
                              start(cor.gr) > beg & 
                              end(cor.gr) < END]
  }
  

  
  if(length(selected_rows) > 0){
    selected_rows_track <- AnnotationTrack(
      cor.gr.subset[selected_rows,],
      name = 'Selected',
      fill = elementMetadata(cor.gr.subset[selected_rows,])$colors,
      #col = '#DCDCDC'
      col = 'transparent'
    )
  } else{
    selected_rows_track <- NULL
  }
  
  
  
  corTrack <- AnnotationTrack(
    cor.gr.subset,
    name = 'Cor',
    fill = elementMetadata(cor.gr.subset)$colors,
    #col = '#DCDCDC'
    col = 'transparent'
  )
  
  if(!is.null(motifs_list)){
    motifs_tracks <- attachMotifs(motifs_list, motifs, chr, beg, END)
    
    motifsTrackList <- lapply(1:length(motifs_tracks), function(x) {
      AnnotationTrack(range = motifs_tracks[[x]],
                      fill = pal_jco()(10)[x],
                      #col = '#DCDCDC',
                      col = 'transparent',
                      name = motifs_list[x])
    })
    
      

  }
  
  if(is.null(motifs_list)){
    plotList <- c(idxTrack, axisTrack, covTrackList, peakTrack, selected_rows_track, corTrack, ENSEMBL_Gviz_GeneTrackRegionObject, NULL)
  } else{
    plotList <- c(idxTrack, axisTrack, covTrackList, peakTrack, selected_rows_track, corTrack, ENSEMBL_Gviz_GeneTrackRegionObject, motifsTrackList)
  }

  
  #plotList <- c(idxTrack, axisTrack, covTrackList)
  
  # highlight <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol",'transcript_start','transcript_end'),
  #                    filters = "refseq_mrna", values = transcriptID, mart= ensembl)
  if(length(transcriptID) > 1){
    highlight <- data.frame(transcript_start = getGtfCoords(gene.symbol,ENSEMBL_hg38_local_fromGTF)[[2]], 
                            transcript_end = getGtfCoords(gene.symbol,ENSEMBL_hg38_local_fromGTF)[[3]])
  } else{
    highlight <- transcript_locations %>% filter(refseq_mrna == transcriptID)
  }
  
  
  print("plotting")
  plotTracks(HighlightTrack(plotList,
                            start = highlight$transcript_start,
                            end = highlight$transcript_end,
                            chr = chr), 
             transcriptAnnotation = "gene",
             #collapseTranscripts = "meta",
             ylim = ylims, col.title = 'black', from = beg, to = END,
             title.width = 3)
  
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
#ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#coverage.files <- paste0("Gviz/", list.files("Gviz"))
color.scheme <- read_csv('Color_Scheme.csv')
#color.scheme$name2 = sapply(color.scheme$name, function(x) paste(strwrap(x,5), collapse="\n"))
coverage.files <- color.scheme$file
strsplits <- str_split(coverage.files, "_")
coverage.names <- sapply(strsplits, function(x) paste(x[c(5,2,3)], collapse = " ")) 
coverage.list <- as.list(coverage.files)
names(coverage.list) <- coverage.names

colors = color.scheme$color

#plotGenomeView(coverage.list = coverage.list)


