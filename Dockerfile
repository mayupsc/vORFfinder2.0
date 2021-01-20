# Version: 0.0.1
FROM ubuntu:18.04
FROM continuumio/miniconda3

MAINTAINER YU MA "yuma@psc.ac.cn"

RUN apt-get install -y libcurl4-openssl-dev libxml2-dev gdebi-core libapparmor1 psmisc supervisor libedit2  libssl-dev libncurses5  r-base r-base-dev wget 

#RUN conda install -c conda-forge wget

RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.7.907-amd64.deb; \
    gdebi -n shiny-server-1.5.7.907-amd64.deb;\
    rm shiny-server-1.5.7.907-amd64.deb

## maybe remove faCount

RUN conda install -c bioconda mafft
RUN conda install -c bioconda blast
RUN conda install -c bioconda phylip
RUN conda install -c conda-forge unzip 

RUN /usr/bin/R -e "install.packages(c('shiny', 'shinyjs', 'shinymaterial','shinycssloaders','ggplot2','digest'), repos='http://mirrors.ustc.edu.cn/CRAN/')"

ADD install_bioC.R /src/install_bioC.R 
RUN /usr/bin/Rscript /src/install_bioC.R 

COPY software/latticeExtra_0.6-28.tar /srv/shiny-server/software/

RUN /usr/bin/R -e "install.packages(c('DT'), repos='http://mirrors.ustc.edu.cn/CRAN/')"
ADD install_bioC2.R /src/install_bioC2.R
RUN /usr/bin/Rscript /src/install_bioC2.R

COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY software/ORFfinder /srv/shiny-server/software/ORFfinder
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

EXPOSE 8383


CMD ["/usr/bin/shiny-server.sh"]



#### !!!!!! attention !!!!!!!! ####
## Do not install R by conda ,it will mask /usr/bin/R path, which will lead to shiny-server error. 
## EXPOSE PORT should be Consistent with PORT in shiny-server.conf and run_docker.sh

