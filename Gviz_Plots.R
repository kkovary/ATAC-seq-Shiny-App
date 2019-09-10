###############################################################################
# VISUALIZATION
###############################################################################
library(Gviz)
library(GenomicRanges)
#library(GenomicFeatures)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

# Get BiomaRt objects
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

plotGenomeView <- function(gene.symbol = GENE, left.extension = LEFT, right.extension = RIGHT, genome = "hg38") {
  biomTrack <- BiomartGeneRegionTrack(genome="hg38", name="ENSEMBL", biomart=ensembl, symbol = gene.symbol)
  axisTrack <- GenomeAxisTrack()
  peakTrack <- 
  plotTracks(list(biomTrack,axisTrack), extend.right = right.extension, extend.left = left.extension, transcriptAnnotation = "name")
}

plotTracks(list(biomTrack,axisTrack), extend.right = 50000, extend.left = 250000, transcriptAnnotation = "name")

knownGene <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene)

knownGene <- UcscTrack(genome="hg38", chromosome = "chr1", from=2000000, to=2100000, name="UCSC Genes", track="knownGene", trackType="GeneRegionTrack",
             rstarts="exonStarts", rends="exonEnds", gene="name", symbol="name",
             transcript="name", strand="strand", fill="#8282d2")

plotTracks(knownGene, chromosome = "chr2",  from=2000000, to=3100000, transcript, geneSymbols = T)
