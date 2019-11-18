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
library(biomaRt)

# colors <- c('#E181F4','#7DCD2C','#F6D7B5','#F9DAFF','#D1E8BA',
#             '#EFE90D','#9404B4','#DC7511','#4C9006','#CB21ED',
#             '#F2A760','#F2B9FF','#B0E57C')

###############################################################################
# Input
###############################################################################

GENE = "SYN1"
SLOP = 250000
GENOME = "hg38"

###############################################################################
# Functions
###############################################################################

getEnsemblCoords <- function(GENE) {
  # Too slow, don't use
  bm <- getBM(attributes = c("transcript_start", "transcript_end"), filters = c("hgnc_symbol"), values = c(GENE), mart = ensembl)
  g.start <- min(bm$transcript_start)
  g.end <- max(bm$transcript_end)
  return(c(g.start,g.end))
}

getGvizCoords <- function(gene.symbol) {
  
  biomTrack <- BiomartGeneRegionTrack(genome="hg38", 
                                      name="ENSEMBL", 
                                      biomart=ensembl, 
                                      symbol = gene.symbol)
  
  g.chr <- Gviz::seqlevels(biomTrack)
  pos <- Gviz::position(biomTrack)
  g.beg <- min(pos)
  g.end <- max(pos)
  
  return(list(g.chr, g.beg, g.end, biomTrack))
}

plotGenomeView <- function(gene.symbol = GENE, slop = SLOP, 
                           genome = "hg38", anno.gr = peaks.gr, 
                           coverage.list, ylims = c(0,100),
                           coords = NULL, chr = NULL,
                           beg = NULL, END = NULL,
                           transcriptID = NULL,
                           cor_cut = NULL,
                           cluster_id = NULL) {
  
  #print("creating track")
  #biomTrack <- BiomartGeneRegionTrack(genome=genome, name="ENSEMBL", biomart=ensembl, symbol = gene.symbol)
  
  print("obtaining coordinates")
  if(is.null(coords)){
    coords <- getGvizCoords(gene.symbol)
  }
  
  if(is.null(beg) & is.null(END)){
    chr <- coords[[1]]
    beg <- coords[[2]]
    END <- coords[[3]]
  }
  
  
  biomTrack <- coords[[4]]
  print(as.character(coords[1:3]))
  
  axisTrack <- GenomeAxisTrack()
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
  
  if(cluster_id == 'All'){
    corTrack <- AnnotationTrack(
      cor.gr[cor.gr$transcript_id == transcriptID & cor.gr$mean.gene.corr >= cor_cut,],
      name = 'Cor',
      fill = 'red',
      col = 'transparent'
    )
  } else{
    corTrack <- AnnotationTrack(
      cor.gr[cor.gr$transcript_id == transcriptID & 
               cor.gr$mean.gene.corr >= cor_cut &
               cor.gr$cluster.name == cluster_id,],
      name = 'Cor',
      fill = 'red',
      col = 'transparent'
    )
  }
  
  
  plotList <- c(axisTrack, covTrackList, peakTrack, corTrack, biomTrack)
  
  print("plotting")
  plotTracks(plotList, transcriptAnnotation = "name", 
             ylim = ylims, col.title = 'black', from = beg, to = END)
  
}

###############################################################################
# MAIN
###############################################################################

# Get BiomaRt objects
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

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


