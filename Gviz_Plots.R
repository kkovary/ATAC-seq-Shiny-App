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

getGvizCoords <- function(TRACK) {
  
  g.chr <- Gviz::seqlevels(TRACK)
  pos <- Gviz::position(TRACK)
  g.start <- min(pos)
  g.end <- max(pos)
  
  return(list(g.chr, g.start, g.end))
}

plotGenomeView <- function(gene.symbol = GENE, slop = SLOP, genome = "hg38", anno.gr = peaks.gr, coverage.list, ylims = c(0,100)) {
  
  print("creating track")
  biomTrack <- BiomartGeneRegionTrack(genome=genome, name="ENSEMBL", biomart=ensembl, symbol = gene.symbol)
  
  print("obtaining coordinates")
  coords <- getGvizCoords(biomTrack)
  chr <- coords[[1]]
  start <- coords[[2]]
  end <- coords[[3]]
  print(as.character(coords))
  
  axisTrack <- GenomeAxisTrack()
  idxTrack <- IdeogramTrack(genome = genome, chromosome = chr)
  
  print("filtering")
  filt.gr <- anno.gr[seqnames(anno.gr)==chr & start(anno.gr) > start - slop & end(anno.gr) < end + slop] 
  
  print("importing peaks")
  peakTrack <- AnnotationTrack(filt.gr, name = "Distal Peaks")
  
  print("importing coverage")
  covTrackList <- lapply(1:length(coverage.list), function(x) {
    DataTrack(range = coverage.list[[x]], genome = genome, 
              type = "histogram", name = color.scheme$name[x], 
              chromosome = chr, start = start, end = end,
              col.histogram = colors[x],
              fill.histogram = colors[x],
              col.axis = "black",
              background.title = colors[x])
  })
  
  
  plotList <- c(axisTrack, covTrackList, peakTrack, biomTrack)
  
  #print("plotting")
  plotTracks(plotList, transcriptAnnotation = "name", 
             ylim = ylims, col.title = 'black')
  
}

###############################################################################
# MAIN
###############################################################################

# Get BiomaRt objects
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#coverage.files <- paste0("Gviz/", list.files("Gviz"))
color.scheme <- read_csv('Color_Scheme.csv')
coverage.files <- color.scheme$file
strsplits <- str_split(coverage.files, "_")
coverage.names <- sapply(strsplits, function(x) paste(x[c(5,2,3)], collapse = " ")) 
coverage.list <- as.list(coverage.files)
names(coverage.list) <- coverage.names

colors = color.scheme$color

#plotGenomeView(coverage.list = coverage.list)


