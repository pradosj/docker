# source("https://bioconductor.org/biocLite.R")

library(Rsamtools)         # biocLite("Rsamtools")
library(GenomicAlignments) # biocLite("GenomicAlignments")
library(GenomicFeatures) 
library(ggplot2)           # install.packages("ggplot2")
library(rtracklayer)       # biocLite("rtracklayer")
library(igraph)
library(visNetwork)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Function to load ASQG graphs created by SGA
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
read.asqg <- function(con) {
  txt <- readLines(con)
  
  vt <- read.table(
    text=txt[grepl("^VT\t",txt)],
    sep="\t",stringsAsFactors=FALSE,
    col.names=c("type","contig","dna","tag")
  )[-1]
  
  ed <- read.table(
    text=txt[grepl("^ED\t",txt)],
    stringsAsFactors=FALSE,
    col.names=c("type","seq1","seq2","start1","end1","len1","start2","end2","len2","dir","overlap_mismatch")
  )[-1]
  
  graph_from_data_frame(ed,vertices=vt)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Convert contigs ranges to genomic ranges given their alignment
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
ctgLocs2refLocs <- function(K,aln) {
  Kcov <- coverage(as(ranges(K),"IRangesList"),width=unname(seqlengths(K)[as.integer(seqnames(K))]))
  Kcov <- RleList(Kcov>0)
  
  h <- findMatches(seqnames(K),names(aln))
  mcols(h)$qcov <- local({
    h.cov <- Kcov[queryHits(h)]
    # TODO: check if this should be after or before remove H/S in cigar
    h.cov[strand(aln[subjectHits(h)])=="-"] <- revElements(h.cov[strand(aln[subjectHits(h)])=="-"])
    h.cov[cigarRangesAlongQuerySpace(cigar(aln[subjectHits(h)]),before.hard.clipping=TRUE,ops=c("H","S"))] <- NA
    h.cov <- h.cov[!is.na(h.cov)]
    IRangesList(h.cov)
  })
  h <- h[lengths(mcols(h)$qcov)>0]
  mcols(h)$qcov <- unlist(mcols(h)$qcov)
  mcols(h)$genomic <- GRanges(seqnames(aln[subjectHits(h)]),pmapFromAlignments(mcols(h)$qcov,aln[subjectHits(h)]),strand(aln)[subjectHits(h)])
  mcols(h)$genomic$cigar_origin <- cigar(aln)[subjectHits(h)]
  mcols(h)$genomic$cigar <- local({
    x <- shift(mcols(h)$genomic,-start(aln)[subjectHits(h)]+1)
    cigarNarrow(mcols(h)$genomic$cigar_origin,start(x),pmin(end(x),cigarWidthAlongReferenceSpace(mcols(h)$genomic$cigar_origin)))
  })
  mcols(h)$genomic$mapq <- mcols(aln)$mapq[subjectHits(h)]
  extractList(mcols(h)$genomic,as(h,"PartitioningByEnd"))
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Compute the variant table from given BAM files
# and annotate with the gff
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
variant.table <- function(bamKmer,bamCtg,gff.file) {
  gff <- import.gff(gff.file,feature.type="gene")
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Load the mapping of the kmers to the contigs and parse informations 
  # contained in headers
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  bam <- readGAlignments(bamKmer, param = ScanBamParam(what=c("qname","seq"),tag = c("NM")))
  mcols(bam)$fwd <- as.integer(sub("^([0-9]*)_([0-9]*)_([0-9]*)$","\\2",mcols(bam)$qname))
  mcols(bam)$rev <- as.integer(sub("^([0-9]*)_([0-9]*)_([0-9]*)$","\\3",mcols(bam)$qname))
  mcols(bam)$kmercount <-  pmin(mcols(bam)$fwd,mcols(bam)$rev)
  aln <- readGAlignments(bamCtg, param = ScanBamParam(what=c("qname","mapq"),tag="NM"),use.names = TRUE)
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Create a GRanges of parts of the contigs covered by a k-mer
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  K <- reduce(unstrand(granges(bam)),with.revmap=TRUE)
  K$numKmer <- lengths(K$revmap)
  K$avgKmerCount <- mean(extractList(mcols(bam)$kmercount,K$revmap))
  K$maxKmerCoverage <- max(RleList(coverage(bam))[K])
  K$revmap <- NULL
  K$coveredGenomicRegions <- ctgLocs2refLocs(K,aln)
  
  K <- local({
    G <- unlist(K$coveredGenomicRegions)
    nn <- nearest(unstrand(G),gff)
    G$nearest_gene_tag <- gff$locus_tag[nn]
    G$nearest_gene_distance[!is.na(nn)] <- distance(unstrand(G),gff[nn[!is.na(nn)]])
    
    K$coveredGenomicRegions <- relist(G,K$coveredGenomicRegions)
    K$coveredGenomicRegionsStr <- unstrsplit(relist(paste0(
      as.character(G),
      ":",G$cigar_origin,
      ":",G$cigar,
      ":",G$mapq,
      ":",G$nearest_gene_tag,
      ":",G$nearest_gene_distance
    ),K$coveredGenomicRegions),";")
    K$rank <- min(relist(rank(G,ties.method="first"),K$coveredGenomicRegions))
    
    K$min_mapq <- min(relist(G$mapq,K$coveredGenomicRegions))
    K$max_cigarI <- max(relist(cigarOpTable(G$cigar)[,"I"],K$coveredGenomicRegions))
    K$max_cigarD <- max(relist(cigarOpTable(G$cigar)[,"D"],K$coveredGenomicRegions))
    K$chr <- unstrsplit(CharacterList(unique(relist(seqnames(G),K$coveredGenomicRegions))),";")
    K$min_chr_start <- min(relist(start(G),K$coveredGenomicRegions))
    K$nearest_gene_tags <- unstrsplit(unique(relist(G$nearest_gene_tag,K$coveredGenomicRegions)),";")
    K$min_nearest_gene_distance <- min(relist(G$nearest_gene_distance,K$coveredGenomicRegions))
    
    K
  })
  K <- K[order(K$rank)]
  K$contig <- as.character(seqnames(K))
  lib <- basename(dirname(bamCtg))
  seqlevels(K) <- paste(lib,seqlevels(K),sep="/")
  K
}



# plot of contigs stats
plot.variants.stats <- function(K) {
  K$within_gene <- ifelse(K$min_nearest_gene_distance==0,"yes","no")
  ggplot(as.data.frame(K),aes(x=maxKmerCoverage,y=avgKmerCount)) + 
    geom_point(aes(size=width,color=chr,shape=within_gene)) + 
    geom_text(aes(label=sub("^merged-","",basename(as.character(seqnames)))),size=3,vjust=2,color="lightgrey") +
    #geom_vline(xintercept=15,linetype=2) +
    #scale_y_log10() + 
    theme_bw()
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load assembly graph
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
visAssemblyGraph <- function(asqg.file,bamCtg,nodes) {
  g <- read.asqg(asqg.file)
  aln <- readGAlignments(bamCtg,use.names = TRUE)
  #G <- induced_subgraph(g,c(V(g)[name%in%nodes | .nei(V(g)[name%in%setdiff(nodes,names(aln))])]))
  G <- induced_subgraph(g,c(V(g)[name%in%nodes | .nei(V(g)[name%in%nodes])]))
  V(G)$shape <- ifelse(V(G)$name %in% nodes,"circle","box")
  
  V(G)$chr <- unstrsplit(unique(CharacterList(split(seqnames(aln),names(aln)))),",")[V(G)$name]
  V(G)$chr[is.na(V(G)$chr)] <- "unmapped"
  V(G)$group <- as.factor(V(G)$chr)
  V(G)$genomic <- unstrsplit(CharacterList(split(as.character(granges(aln)),names(aln))),"\n")[V(G)$name]
  V(G)$name <- paste0(V(G)$name,"\n",nchar(V(G)$dna),"bp\n",V(G)$genomic)
  
  visIgraph(G) %>%
    visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>%
    visLegend()
} 



# excel table
export.variant.tsv <- function(K,tsv.file) {
  K$num_covered_genomic_region <- lengths(K$coveredGenomicRegions)
  K$contig_len <- seqlengths(K)[as.character(seqnames(K))]
  K <- as.data.frame(K)
  cn <- c("seqnames","contig","contig_len","start","end","width",
          "numKmer","avgKmerCount","maxKmerCoverage",
          "num_covered_genomic_region",
          "min_mapq","max_cigarI","max_cigarD","chr","min_chr_start",
          "nearest_gene_tags","min_nearest_gene_distance",
          "coveredGenomicRegionsStr")
  K <- K[intersect(cn,colnames(K))]
  write.table(K,file=tsv.file,row.names=FALSE,sep="\t")
}

read.variant.table <- function(var.file) {
  read.table(var.file,sep="\t",header=TRUE,stringsAsFactors=FALSE)  
}






