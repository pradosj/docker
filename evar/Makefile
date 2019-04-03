
-include Makefile.config

REF ?= /tmp/SA564.fa
REF_GFF ?= /tmp/SA564.gff
SGA_PREPROC_FLAGS ?= -r AGATCGGAAGAGCACACGTCTGAACTCCAGTC -c AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
SGA_FILTER_FLAGS ?= -x10 -k35
SGA_FMMERGE_FLAGS ?= -m35
SGA_KMERCOUNT_FLAGS ?= -k35
SGA_THREAD ?= -t4
SGA_OVERLAP_FLAGS ?= -m35
KMER_FILTER_THRESHOLD ?= 10
BWA_KMERMAP_FLAGS ?= -T10
BACKGROUND_BWT ?= $(REF:%.fa=%).bwt

.PRECIOUS: %.bwamem.bam %.bwt %.fa.pac %.preproc.bwt %.filter.pass.fa %.preproc.fa %.rmdup.merged.fa


%.all:
	$(MAKE) $*.preproc.filter.pass.rmdup.merged.bwamem.bam \
					$*.preproc.filter.pass.rmdup.merged.bwamem.bam.bai \
	        $*.bwamem.bam $*.bwamem.bam.bai \
	        $*.preproc.filter.pass.kmercount.tsv \
	        $*.preproc.filter.pass.kmercount.fa \
	        $*.preproc.filter.pass.rmdup.merged.fa \
	        $*.preproc.filter.pass.rmdup.merged.kmermap.bam \
	        $*.preproc.filter.pass.kmercount.bwamem.bam \
	        $*.preproc.filter.pass.kmercount.bwamem.bam.bai \
	        $*.preproc.filter.pass.rmdup.merged.asqg.gz \
	        $*.preproc.filter.pass.rmdup.merged.variants.tsv \
	        $*.preproc.filter.pass.rmdup.merged.variants.bed \
	        $*.preproc.filter.pass.rmdup.merged.variants.pdf \
	        $*.preproc.filter.pass.rmdup.merged.variants.html
 
%.preproc.bwt : %.preproc.fa
	cd $(@D);sga index $(SGA_THREAD) -a ropebwt $(notdir $^)

%.bwt : %.fa
	cd $(@D);sga index $(SGA_THREAD) $(notdir $^)

%.fa.pac : %.fa
	bwa index $^

%.bwamem.bam : $(REF).pac %_R1.fastq.gz %_R2.fastq.gz
	bwa mem $(SGA_THREAD) $(REF) $(word 2,$^) $(word 3,$^) | samtools view -Sb - | samtools sort -o $@ -

%.bwamem.bam : $(REF).pac %.fastq.gz
	bwa mem $(SGA_THREAD) $(REF) $(word 2,$^) | samtools view -Sb - | samtools sort -o $@ -

%.bwamem.bam : $(REF).pac %.fa 
	bwa mem $(SGA_THREAD) $(REF) $(word 2,$^) | samtools view -Sb - | samtools sort -o $@ -

%.bam.bai : %.bam
	samtools index $^

%.preproc.fa : %_R1.fastq.gz %_R2.fastq.gz
	cd $(@D);sga preprocess -p1 -q3 -f1 $(SGA_PREPROC_FLAGS) $(notdir $^) > $(notdir $@)

%.preproc.fa : %.fastq.gz
	cd $(@D);sga preprocess -p0 -q3 -f1 $(SGA_PREPROC_FLAGS) $(notdir $^) > $(notdir $@)

%.filter.pass.fa %.filter.pass.bwt: %.fa %.bwt
	cd $(@D);sga filter $(SGA_THREAD) $(SGA_FILTER_FLAGS) --kmer-both-strand --no-duplicate-check $(notdir $<)

%.rmdup.fa : %.fa
	cd $(@D);sga rmdup $(SGA_THREAD) $(notdir $^)

%.rmdup.merged.fa : %.rmdup.fa %.rmdup.bwt
	cd $(@D);sga fm-merge $(SGA_THREAD) $(SGA_FMMERGE_FLAGS) $(notdir $<)

%.kmercount.tsv : %.bwt $(BACKGROUND_BWT)
	sga kmer-count $(SGA_KMERCOUNT_FLAGS) $^ |awk '($$2>$(KMER_FILTER_THRESHOLD)) && ($$3>$(KMER_FILTER_THRESHOLD)) && ($$4==0) && ($$5==0)' > $@

%.kmercount.fa : %.kmercount.tsv
	cat $^ | awk '{print ">" NR "_" $$2 "_" $$3; print $$1}' > $@

%.rmdup.merged.kmermap.bam : %.rmdup.merged.fa.pac %.kmercount.fa 
	bwa mem -a $(BWA_KMERMAP_FLAGS) $*.rmdup.merged.fa $*.kmercount.fa | awk '($$6~/^[0-9]+M$$/ && $$12~/NM:i:0/) || $$0~/^@/' |samtools view -Sb - | samtools sort -o $@ -

%.asqg.gz : %.fa %.bwt
	cd $(@D);sga overlap $(SGA_THREAD) $(SGA_OVERLAP_FLAGS) $(notdir $<) 

%.clean :
	cd $(@D);rm -f $**.bwt $**.rbwt $**.sai $**.rsai $**.sa $**.pac $**.ann $**.amb $*.preproc.fa $*.preproc.discard.fa $*.preproc.filter.pass.fa $*.preproc.filter.pass.rmdup.dups.fa $*.preproc.filter.pass.rmdup.fa $(REF).rsai $(REF).rbwt $(REF).sai $(REF).sa $(REF).amb $(REF).pac $(REF).bwt $(REF).ann $(REF:%.fa=%).sa $(REF:%.fa=%).pac $(REF:%.fa=%).bwt $(REF:%.fa=%).rsai $(REF:%.fa=%).rbwt $(REF:%.fa=%).sai

%.variants.bed %.variants.tsv : %.bwamem.bam %.kmermap.bam $(REF_GFF)
	Rscript -e 'source("/tmp/Analyse_bam.R");vartbl <- variant.table("$(word 2,$^)","$(word 1,$^)","$(REF_GFF)");export.variant.tsv(vartbl,"$*.variants.tsv");export.bed(reduce(unlist(vartbl$$coveredGenomicRegions)),"$*.variants.bed")'

%.variants.pdf : %.variants.tsv
	Rscript -e 'source("/tmp/Analyse_bam.R");vartbl <- read.variant.table("$<");ggsave("$@",plot.variants.stats(vartbl))'

%.variants.html : %.variants.tsv %.asqg.gz %.bwamem.bam
	Rscript -e 'source("/tmp/Analyse_bam.R");vartbl <- read.variant.table("$<");visSave(visAssemblyGraph("$(word 2,$^)","$(word 3,$^)",vartbl[["contig"]]),"$@")'

