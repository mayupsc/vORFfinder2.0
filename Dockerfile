# Version: 0.0.1
FROM ubuntu:18.04
FROM continuumio/miniconda3

MAINTAINER YU MA "yuma@psc.ac.cn"

RUN apt-get install -y libcurl4-openssl-dev libxml2-dev gdebi-core libapparmor1 psmisc supervisor libedit2  libssl-dev libncurses5

RUN conda install -c r r-base=3.6.1
RUN conda install -c conda-forge wget

## maybe remove faCount
RUN conda install -c bioconda ucsc-facount
RUN conda install -c bioconda mafft
RUN conda install -c bioconda blast
RUN conda install -c bioconda phylip
RUN conda install -c r r-ggplot2
RUN conda install -c conda-forge unzip 



RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.7.907-amd64.deb; \
    gdebi -n shiny-server-1.5.7.907-amd64.deb;\
    rm shiny-server-1.5.7.907-amd64.deb


RUN R -e "install.packages(c('shiny', 'shinyjs', 'shinymaterial'), repos='http://mirrors.ustc.edu.cn/CRAN/')"

ADD install_bioC.R /src/install_bioC.R 

RUN Rscript /src/install_bioC.R 


EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY software/ORFfinder /srv/shiny-server
CMD ["/usr/bin/shiny-server.sh"]
