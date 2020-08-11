FROM bioconductor/release_core2

RUN apt-get update && apt-get install -y pigz
RUN Rscript -e 'BiocManager::install(c("Rbowtie","ShortRead"),suppressUpdates="")'
RUN Rscript -e 'install.packages("igraph")'
#RUN gpg --keyserver keyring.debian.org --recv-keys 04EE7237B7D453EC 648ACFD622F3D138 

ADD Makefile src /tmp/
VOLUME /export/
WORKDIR /export/
