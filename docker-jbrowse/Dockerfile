FROM httpd

RUN apt-get update && apt-get install -y \     
    unzip \
    build-essential \
    libpng-dev \
    zlib1g-dev \
    libgd2-xpm-dev 

COPY JBrowse-1.12.0.zip /tmp/
RUN unzip -d /usr/local/apache2/ /tmp/JBrowse-1.12.0.zip && rm -rf /tmp/JBrowse-1.12.0.zip /usr/local/apache2/htdocs && mv /usr/local/apache2/JBrowse-1.12.0 /usr/local/apache2/htdocs
WORKDIR /usr/local/apache2/htdocs/
RUN ./setup.sh


