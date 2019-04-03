# evar


Variants calling pipieline by de novo assembly using SGA, BWA, Samtools and R script for pool sequencing

---------------------------------------------------------------------------------------------------------------------------------

files needed

	reference genome in .fa 
	genomic annotation in .gff
	single-end one file fastq.gz or double-end two files R1.fastq.gz R2.fastq.gz

	! your reference can't have any N in place of A, T, G or C --> sript will crash

---------------------------------------------------------------------------------------------------------------------------------

default parameters

	# -r = remove-adapter-fwd (str)
	# -c = remove-adapter-rev (str)
	SGA_PREPROC_FLAGS:= -r AGATCGGAAGAGCACACGTCTGAACTCCAGTC -c AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG

	# -x = min x kmer coverage in a read (int)
	# -k = size of kmer (int)
	SGA_FILTER_FLAGS:= -x15 -k61

	# -m = minimum overlap required between two reads to merge (int)
	SGA_FMMERGE_FLAGS:= -m61

	# -k = size of kmer (int)
	SGA_KMERCOUNT_FLAGS:= -k61

	# -t = number of thread (int)
	SGA_THREAD:= -t4

	# -m = minimum overlap required between two reads to merge (int)
	SGA_OVERLAP_FLAGS:= -m61

  Create Makefile.config in your working directory with the following template
  
  ```
  REF = ref.fa
  REF_GFF = ref.gff

  SGA_PREPROC_FLAGS:= -r AGATCGGAAGAGCACACGTCTGAACTCCAGTC -c AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
  SGA_FILTER_FLAGS:= -x15 -k61
  SGA_FMMERGE_FLAGS:= -m61
  SGA_KMERCOUNT_FLAGS:= -k61
  SGA_THREAD:= -t4
  SGA_OVERLAP_FLAGS:= -m61
  ```

	Some parameters of rules were not in parameters but can be change in the original makefile
	For more details on parameters see the help for SGA, BWA and samtools

---------------------------------------------------------------------------------------------------------------------------------

docker command to run the application

	```
	docker run --rm -v $(pwd):/export benjamindn/evar-call fastq.all
	```

---------------------------------------------------------------------------------------------------------------------------------

IGV visualisation

	#all your reads are mapped
	.bwamem.bam

	#De novo assembly 
	.preproc.filter.pass.kmercount.bwamem.bam

	#Region of interest by kmer
	.preproc.filter.pass.rmdup.merged.bwamem.bam

	The best way to visualize your result is to put together the de novo assembly and the region of interest in IGV

---------------------------------------------------------------------------------------------------------------------------------

Example

	Create a directory data 
	Put your fastq.gz, reference.fa and annotation.gff in this directory
	Go in the data directory
	Run this command 
	```
	docker run --rm -v $(pwd):/export benjamindn/evar-call fastq.all
	```
	Wait for all the processus (more your genome to analyse is big more times is needed)
	At the end you get the three files for an IGV visualisation, a file in .tsv you can open with a tablesoftware, a repartition plot of mutation and an interactive assembly network to see the de novo assembly.  

---------------------------------------------------------------------------------------------------------------------------------

Table explanation

	work in progress
