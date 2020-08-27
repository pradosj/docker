

The current version include:

* samtools 1.3.1
* bcftools 1.3.1
* htslib 1.3.2
* bwa 0.7.15
* bowtie 2.2.9
* tophat 2.1.1
* cuffdiff 2.1.1
* sga 0.10.15
* STAR 2.5.2b
* picard-tools
* HTSeq 
* MACS2




Example usage to map contigs onto a reference genome with `BWA mem`.
```
docker run --rm -v $(pwd):/export pradosj/ngs 'bwa index ref.fa && bwa mem -t4 ref.fa contigs.fasta | samtools view -Sb - | samtools sort -o contigs.bam - && samtools index contigs.bam'
```






