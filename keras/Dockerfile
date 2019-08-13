FROM bioconductor/release_core2

RUN Rscript --vanilla -e 'BiocManager::install("rhdf5")'
RUN Rscript --vanilla -e 'BiocManager::install("HDF5Array")'
RUN apt-get update && apt-get install -y \
    libssl-dev \
    python3 \
    python3-pip \
    python3-virtualenv && \
    rm -rf /var/lib/apt/lists/*
RUN Rscript --vanilla -e 'install.packages(c("keras","ggplot2"),repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'keras::install_keras()'


