FROM ubuntu

RUN apt-get update && apt-get install -y \
    default-jdk \
    unzip

ADD http://data.broadinstitute.org/igv/projects/downloads/IGV_2.3.88.zip /
WORKDIR /
RUN unzip IGV_2.3.88.zip

VOLUME ["/export/"]
ENTRYPOINT ["java","-Xmx4000m","-Dapple.laf.useScreenMenuBar=true","-Djava.net.preferIPv4Stack=true","-jar","/IGV_2.3.88/igv.jar"]

