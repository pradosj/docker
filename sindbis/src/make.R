


#-#-#-#-#-#-#-#-#-#
# WHICH REGION IS EACH READ COMING FROM?
# build the index mapping file that associates to any sequence an index (corresponding to a region)
#-#-#-#-#-#-#-#-#-#
make.indexmap <- function(index.file,indexmap.file) {
  library(IRanges)
  library(Biostrings)
  
  index <- as(read.table(index.file,sep="\t",header=TRUE,stringsAsFactors = FALSE),"DataFrame")
  index$sequence <- gsub(" *","",index$sequence)
  index$sequence <- reverseComplement(DNAStringSet(index$sequence))
  
  map <- DataFrame(seq=DNAStringSet(mkAllStrings(c("A","C","G","T","N"),6)))
  n <- sapply(DNAStringSet(index$sequence),neditAt,map$seq,at=1)
  map$closest.index.id <- max.col(-n)
  map$closest.index.seq <- index$sequence[map$closest.index.id]
  map$closest.index.name <- index$name[map$closest.index.id]
  map$closest.index.mismatch <- n[cbind(seq_along(map$closest.index.id),map$closest.index.id)]
  
  map
}


#-#-#-#-#-#-#-#-#-#
# map a FASTA file onto the given bowtie index (Justus library sequences???), and generate a BAM
# The method uses Bowtie and find all possible alignments with 3 mismatches allowed
#-#-#-#-#-#-#-#-#-#
bowtie.map <- function(fa.file,bwt.index,bam.file) {
  library(Rbowtie)
  library(Rsamtools)
  
  tmp.sam <- tempfile(fileext=".sam")
  bowtie(sequences=fa.file,index=bwt.index,outfile=tmp.sam,f=TRUE,p=3,n=3,v=3,l=32,sam=TRUE,all=TRUE,best=TRUE,force=TRUE) # v = allowed mismatches; n = max mismatches in seed (0 to 3, if seed = 28 nt); p = number of CPU for analysis (keep 3) ; l = length of the seed (28 by default)
  asBam(tmp.sam,sub(".bam$","",bam.file),overwrite=TRUE)
  return(bam.file)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Find clusters in a self-mapped BAM graph
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
make.barcode.cluster <- function(self.bam) {
  library(igraph)
  library(GenomicAlignments)
  library(Biostrings)
  library(S4Vectors)
  
  # Load the barcode graph
  aln <- readGAlignments(self.bam,param=ScanBamParam(what="qname",tag="NM",flag = scanBamFlag(isMinusStrand = FALSE)))
  g <- graph_from_edgelist(cbind(mcols(aln)$qname,as.character(seqnames(aln))),directed=FALSE)
  #E(g)$weight <- mcols(aln)$NM
  
  #-#-#-#-#-#-#-#-#-#
  # Cluster the barcodes
  # the cluster center is the most abundant sequence within a cluster, considered as the "real" sequence of this barcode
  # 3 mismatches are allowed between barcodes to define them as bellonging to a same cluster
  #-#-#-#-#-#-#-#-#-#
  k <- components(g)
  k <- DataFrame(vertex=names(k$membership),cluster.id=unname(k$membership),cluster.size=k$csize[k$membership])
  k$barcode.count <- as.integer(sub("[A-Z]*_","",k$vertex))
  k$barcode <- DNAStringSet(sub("_[0-9]*","",k$vertex))
  o <- order(k$barcode.count,decreasing=TRUE)
  k$is.cluster.center[o] <- !duplicated(k$cluster.id[o])
  k$cluster.barcode.count <- as.vector(tapply(k$barcode.count,k$cluster.id,sum)[k$cluster.id])
  
  k
}


#-#-#-#-#-#-#-#-#-#
# Map all barcodes of a pup onto filtered S1 clusters identified
#-#-#-#-#-#-#-#-#-#
make.cluster.assignment <- function(in.tsv,in.bam,out.tsv.gz) {
  library(GenomicAlignments)
  aln <- readGAlignments(in.bam,param=ScanBamParam(what="qname",tag="NM",flag=scanBamFlag(isMinusStrand=FALSE)))
  x <- read.table(in.tsv,sep="\t",stringsAsFactors=FALSE,header=TRUE)
  i <- match(x$bc32,mcols(aln)$qname)
  x$barcode.mapping <- as.character(seqnames(aln))[i]
  x$barcode.mapping.err <- mcols(aln)$NM[i]
  con <- gzfile(out.tsv.gz,"w")
  on.exit(close(con))
  write.table(x,file=con,sep="\t",quote=FALSE,row.names=FALSE,na="")
}





