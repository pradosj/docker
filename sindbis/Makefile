
SHELL=bash


#test command:
#docker run --rm -v $(pwd):/export sindbis -j2 MOS5_P7_{S1,RN}.umi.class.viral.fasta
#docker run --rm -v $(pwd):/export sindbis MOS5_P7_S1.umi.class.viral{.bowtie_index/index.1.ebwt,.self.bam.clusters.filter100.fasta}
#docker run --rm -v $(pwd):/export sindbis MOS5_P7_S1.umi.class.viral.self.bam.clusters.filter100{.bowtie_index/index.1.ebwt,.bowtie_aln/MOS5_P7_RN.umi.class.viral.bam}


.PHONY:usage
usage:
	@echo 'usage: docker run --rm -v $$(pwd):/export pradosj/sindbis <file.sindbis/all>'
	@echo '  <file.fastq.gz> and <file.index.tsv> must exist'
	@echo 'Run reads processing pipeline for a MAPseq experiment'


.PHONY:%.sindbis/all
%.sindbis/all:%.sindbis/demux_report.txt %.sindbis/links.txt
	awk '{ \
	  ref = $$1 "_" $$2 ".umi.class.viral.self.bam.clusters.filter100"; \
	  targ = $$1 "_" $$3 ".umi.class.viral"; \
	  print "$(@D)/" ref ".fasta"; \
	  print "$(@D)/" ref ".bowtie_index/index.1.ebwt"; \
	  print "$(@D)/" targ ".fasta"; \
	  print "$(@D)/" ref ".bowtie_aln/" targ ".bam"; \
	  print "$(@D)/" ref ".bowtie_aln/" targ ".bam.tsv.gz"; \
	}' $*.sindbis/links.txt | xargs $(MAKE) -f /tmp/Makefile 



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Demultiplexing & sequence parsing
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# optionally cut reads to a size of 105bp
%.105bp.fastq.gz:%.fastq.gz
	gzip -dc $< | awk 'NR%2==0{$$0=substr($$0,1,105)}{print}' | pigz > $@

# Generate demultiplexing table
# Reverse complement the index sequence and only keep perfectly matching indexes
.PRECIOUS:%.sindbis/demuxmap.tsv
%.sindbis/demuxmap.tsv:%.index.tsv
	mkdir -p $(@D) && Rscript -e 'source("/tmp/src/make.R");map <- make_demuxmap("$<");map <- subset(map,closest_index_mismatch<=0);write.table(map[c("seq","closest_index_name")],"$@",sep="\t",quote=FALSE,row.names=FALSE)'

# Demultiplex FASTQ reads according to demuxmap
.PRECIOUS:%.sindbis/demux_report.txt %.sindbis/dissections.txt %.sindbis/links.txt
%.sindbis/demux_report.txt:%.sindbis/demuxmap.tsv %.fastq.gz
	gzip -dc $(word 2,$^) | awk -v demuxmap='$(word 1,$^)' -v odir='$(@D)' -f /tmp/src/fq_demux.awk 

# Create dissections.txt with a line for each dissection
%.sindbis/dissections.txt:%.sindbis/demux_report.txt
	awk '$$1~/.*_.*.fastq.gz/{ \
	  sub(".fastq.gz$$","",$$1); \
	  targ = pup = $$1; \
	  sub("_.*","",pup); \
	  targ = substr(targ,length(pup)+2); \
	  print pup "\t" targ; \
	}' $< | sort -k1,1 -k2,2 > $@

# Create links.txt discribing links between dissections
%.sindbis/links.txt:%.sindbis/dissections.txt
	join <(awk '$$2~/Inj/' $<) <(awk '$$2!~/Inj/' $<) > $@


# Generate UMI depuplicated table of the Barcodes found in the FASTQ reads
# ignore reads with an N in the UMI or in the barcode
%.umi.tsv.gz:%.fastq.gz
	gzip -dc $(word 1,$^) | \
	awk 'BEGIN{OFS="\t"}(NR%4==2){umi=substr($$0,95,12);bc32=substr($$0,1,32);if ((umi!~/N/) && (bc32!~/N/)){print umi,bc32;}}' | \
	LC_ALL=C sort --parallel=6 --buffer-size=75% --compress-program=pigz | uniq -c | \
	awk 'BEGIN{OFS="\t";print "umi","bc32","umi_freq";}{print $$2,$$3,$$1}' | pigz > $@
	

# Classify reads accoring to the content of the bc32 barcode
# Allow 1 mismatch in the spike sequence
.PRECIOUS:%.umi.class.tsv.gz
%.umi.class.tsv.gz:%.umi.tsv.gz
	gzip -dc $(word 1,$^) | \
	awk 'BEGIN{ \
	       OFS="\t"; \
	       print "umi","bc32","umi_freq","class"; \
	       yy["CC"]=yy["TT"]=yy["CT"]=yy["TC"]=0; \
	       spike["AACAGTCA"]=spike["ACCAGTCA"]=spike["AGCAGTCA"]=spike["ATAAGTCA"]=spike["ATCAATCA"]=spike["ATCACTCA"]=spike["ATCAGACA"]=spike["ATCAGCCA"]=spike["ATCAGGCA"]=spike["ATCAGTAA"]=spike["ATCAGTCA"]=spike["ATCAGTCC"]=spike["ATCAGTCG"]=spike["ATCAGTCT"]=spike["ATCAGTGA"]=spike["ATCAGTTA"]=spike["ATCATTCA"]=spike["ATCCGTCA"]=spike["ATCGGTCA"]=spike["ATCTGTCA"]=spike["ATGAGTCA"]=spike["ATTAGTCA"]=spike["CTCAGTCA"]=spike["GTCAGTCA"]=spike["TTCAGTCA"]=0; \
	     } \
	     (NR>1){ \
	       cl="unknown";bc32=$$2; \
	       if (substr(bc32,25,8) in spike){cl="spike"} else {if (substr(bc32,31,2) in yy) cl="viral"}; \
	       print $$1,$$2,$$3,cl; \
	     }' | pigz > $@

.PRECIOUS:%.umi.class.tsv.gz.stat
%.umi.class.tsv.gz.stat:%.umi.class.tsv.gz
	gzip -dc $< | awk '($$4!="class"){n[$$4]++}END{for(i in n) print i "\t" n[i]}' > $@

# extract a FASTA with unique viral barcodes from a umi.class.tsv.gz file
.PRECIOUS:%.umi.class.viral.fasta
%.umi.class.viral.fasta:%.umi.class.tsv.gz
	gzip -dc $(word 1,$^) | awk '$$NF=="viral"{print $$2}' | sort | uniq -c | awk '$$1>=1{print ">" $$2 "_" $$1 "\n" $$2}' > $@



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Injection-sites self mapping for barcode error correction
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# index a FASTA sequence with bowtie_index
.PRECIOUS:%.bowtie_index/index.1.ebwt
%.bowtie_index/index.1.ebwt:%.fasta
	 Rscript -e 'library("Rbowtie");bowtie_build("$<","$(@D)")'

# Map a fasta against himself
.PRECIOUS:%.self.bam
%.self.bam:%.fasta
	$(MAKE) -f /tmp/Makefile $*.bowtie_index/index.1.ebwt
	Rscript -e 'source("/tmp/src/make.R");bowtie_map("$(word 1,$^)","$*.bowtie_index/index","$@")'

# Find sequence clusters from a BAM generated by self-mapping a FASTA 
.PRECIOUS:%.self.bam.clusters.tsv
%.self.bam.clusters.tsv:%.self.bam
	Rscript -e 'source("/tmp/src/make.R");k <- make_barcode_cluster("$<");write.table(k,"$@",sep="\t",quote=FALSE,row.names=FALSE)'

# Generate a FASTA with the centers of the clusters that have >=100 barcodes
.PRECIOUS:%.self.bam.clusters.filter100.fasta
%.self.bam.clusters.filter100.fasta:%.self.bam.clusters.tsv
	awk 'BEGIN{FS="\t"} ($$7>=100 && $$6=="TRUE"){print ">" $$1 "_" $$3 "_" $$7;print $$5}' $< > $@





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Barcode mapping onto a reference library of barcodes
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Map any fasta onto a given ref
.PRECIOUS:%.bam
%.bam:
	$(MAKE) $(@D:%.bowtie_aln=%.bowtie_index)/index.1.ebwt 
	mkdir -p $(@D) && $(MAKE) $(@D)/../$(notdir $*).fasta
	Rscript -e 'source("/tmp/src/make.R");bowtie_map("$(@D)/../$(notdir $*).fasta","$(@D:%.bowtie_aln=%.bowtie_index)/index","$@")'

%.viral.bam.tsv.gz:%.viral.bam
	$(MAKE) $(@D)/../$(notdir $*).tsv.gz
	Rscript -e 'source("/tmp/src/make.R");make_cluster_assignment("$(@D)/../$(notdir $*).tsv.gz","$(word 1,$^)","$@")'




