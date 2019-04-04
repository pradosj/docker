

EMOTE is a biological protocol that rely on next generation sequencing to detect 5' ends of transcripts on a global scale. This docker image provides necessary software to analyze sequencing data generated with EMOTE protocol. The image inherits bioconductor/release_core which contains an R installation, core Bioconductor packages, and Rstudio server for the UI. Additionaly, it installs the following R packages:

* roxygen2
* VGAM
* ShortRead
* Rbowtie
* EMOTE

To run the image use:
```
docker run -d -p 8787:8787 -v $(pwd):/export \
           -e USER=$USER -e PASSWORD=EMOTE \
           pradosj/docker-emote
```


