FROM bioconductor/release_core2

RUN Rscript --vanilla -e 'BiocManager::install("GenomicFiles")'

RUN apt-get update && apt-get install -y \
    curl \
    libncurses-dev
    
#-#-#-#-#-#-#-#-#-#-#-#-#
# SAMTOOLS
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/htslib-1.3.2 && make && make install && \
    rm -rf /tmp/htslib-1.3.2
RUN curl -kL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/samtools-1.3.1 && make && make install && \
    rm -rf /tmp/samtools-1.3.1
RUN curl -kL https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bcftools-1.3.1 && make && make install && \
    rm -rf /tmp/bcftools-1.3.1

#-#-#-#-#-#-#-#-#-#-#-#-#
# BWA
#-#-#-#-#-#-#-#-#-#-#-#-#
RUN curl -kL http://netix.dl.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2 | tar -C /tmp -jxf - && \
    cd /tmp/bwa-0.7.15 && make && find /tmp/bwa-0.7.15/ -type f -executable -exec mv '{}' /usr/local/bin/ ';' && \
    rm -rf /tmp/bwa-0.7.15/




