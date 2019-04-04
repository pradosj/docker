# evar

**evar** is a variants calling pipieline by *de novo* assembly using SGA, BWA, Samtools and R scripts.


### Required input files

The pipeline require several input files:

* sequenced reads in FASTQ format:
    * either single-end reads with file extension `.fastq.gz`
    * either paired-end reads with file extension `_R1.fastq.gz` and `_R2.fastq.gz`
* reference genome in FASTA format, with file extension `.fa`. This reference sequence can't contains N but only A,T,G,C
* genomic annotation in GFF format, with extension `.gff`
	


### Quick start

1. Install **docker** from https://docker.com

2. Create a new directory:
```
mkdir my_project
cd my_project
```

3. Put your reference genome (FASTA format) and genes annotations (GFF format) inside: 
```
cp ref.fa ref.gff ./
```
**WARNING:** No N are accepted in the FASTA file otherwise the pipeline crashes

4. Put your sequenced reads inside, in gzipped-FASTQ format:
```
cp reads_R1.fastq.gz reads_R2.fastq.gz ./
```

5. Call the variants
```
docker run --rm -v $(pwd):/export plinderlab/evar KMER=61 KMER_NOISE_LEVEL=10 REF_FA=ref.fa REF_GFF=ref.gff reads.varcall
```
The process is quite time consuming so be patient !


### Output files

IGV visualisation

	#all your reads are mapped
	.bwamem.bam

	#De novo assembly 
	.preproc.filter.pass.kmercount.bwamem.bam

	#Region of interest by kmer
	.preproc.filter.pass.rmdup.merged.bwamem.bam


	The best way to visualize your result is to put together the de novo assembly and the region of interest in IGV

	At the end you get the three files for an IGV visualisation, a file in .tsv you can open with a tablesoftware, a repartition plot of mutation and an interactive assembly network to see the de novo assembly.  




### Running the variant caller

1. Generate the Makefile
```
docker run --rm -v $(pwd):/export plinderlab/evar
```
2. Edit parameters in Makefile header with a text editor
3. Call the variants
```
docker run --rm -v $(pwd):/export plinderlab/evar reads.varcall
```





