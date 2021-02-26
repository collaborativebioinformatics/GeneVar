FROM ubuntu:18.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true
ARG THREADS=4

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        wget \
        curl \
        less \
        gcc \ 
        samtools \
        tabix \
        tzdata \
        make \
        bzip2 \
        git \
        sudo \
        pkg-config \
	libxml2-dev \
	libssl-dev \
	libmariadbclient-dev \
	libcurl4-openssl-dev \ 
	apt-transport-https \
	software-properties-common \
	dirmngr \
	gpg-agent \ 
        && rm -rf /var/lib/apt/lists/*

ENV TZ=America/Los_Angeles

WORKDIR /build

## R and pandoc
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        && apt-get update \
        && apt-get install -y r-base r-base-dev

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz')"

RUN R -e "install.packages(c('dplyr', 'knitr', 'ggplot2', 'rmarkdown', 'DT', 'ggrepel', 'tidyr', 'BiocManager', 'plotly', 'DiagrammeR', 'visNetwork', 'shiny'))"

RUN R -e "BiocManager::install(c('VariantAnnotation'))"

## more packages
RUN R -e "install.packages(c('rsconnect', 'shinydashboard', 'stringdist'))"

RUN R -e "BiocManager::install(c('Rsamtools'))"
