FROM bioconductor/bioconductor_docker:RELEASE_3_10

RUN apt-get update && apt-get install -y pigz
RUN Rscript -e 'install.packages(c("igraph","devtools","abind","gplots"))'
RUN Rscript -e 'BiocManager::install(c("SummarizedExperiment","rhdf5","HDF5Array"),suppressUpdates="")'

# jpeg & rjson (dependencies for ?)
RUN apt-get update && apt-get install -y libjpeg-dev
RUN Rscript -e 'install.packages(c("jpeg","rjson"))'

# tools 
RUN Rscript -e 'install.packages(c("FNN", "princurve","irlba"))'

# Keras
RUN pip install --upgrade virtualenv
RUN Rscript -e 'install.packages("keras")'
RUN sudo -u rstudio Rscript -e 'keras::install_keras()'

# BMRM
RUN Rscript -e 'devtools::install_github("pradosj/bmrm/bmrm")'

# UMAP
RUN sudo -u rstudio Rscript -e 'reticulate::virtualenv_install(packages="umap-learn")'

# Seurat
RUN Rscript -e 'BiocManager::install("multtest");install.packages("Seurat")'

# Velocyto
RUN Rscript -e 'BiocManager::install("pcaMethods")'
RUN apt-get update && apt-get install -y libhdf5-dev libboost-system-dev libboost-filesystem-dev
RUN Rscript -e 'devtools::install_github("velocyto-team/velocyto.R")'




# More packages
RUN Rscript -e 'install.packages("rmarkdown")'
RUN Rscript -e 'BiocManager::install("rtracklayer")'

