FROM bioconductor/release_core2

RUN Rscript --vanilla -e 'BiocManager::install("rhdf5")'
RUN Rscript --vanilla -e 'BiocManager::install("HDF5Array")'
RUN apt-get update && apt-get install -y \
    libssl-dev \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*
RUN pip3 install tensorflow
RUN pip3 install SciPy
ENV RETICULATE_PYTHON /usr/bin/python3
RUN Rscript --vanilla -e 'install.packages("tensorflow",repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("keras",repos="https://cloud.r-project.org")'

