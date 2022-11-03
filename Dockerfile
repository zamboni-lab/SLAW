FROM rocker/r-ubuntu:22.04

# We consider the values of the data eventually
# The in house R packages are stored in the in_house repository
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN add-apt-repository universe
RUN apt-get install -y build-essential \
  libcurl4-openssl-dev \
  libxml2-dev \
  libssl-dev \
  libzmq3-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libfontconfig1-dev \
  openjdk-8-jre \
  libnetcdf-dev \
  libgmp3-dev \
  software-properties-common \
  && rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# Installing all the dependencies in 1 layer for efficiency
# IF MISSING add ebase and rbase-dev
RUN apt-get -y --no-install-recommends install --fix-missing python3.9 python3-dev gcc gfortran musl-dev g++ python3-pip

RUN apt-get -y --no-install-recommends --fix-missing install openms

# Fix for singularity on HPC
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

# Python libraries
# Always installing the wheel pip3 package first
RUN pip3 install --no-cache-dir wheel psutil
RUN pip3 install numpy
RUN pip3 install --no-cache-dir panda PyYAML
RUN pip3 install sklearn statsmodels

# R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update && apt-get -y --no-install-recommends install --fix-missing r-base r-base-dev
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/"
RUN R -e "install.packages('devtools')"

#Dependencies copy
COPY MZmine-2.52-Linux /MZmine-2.52-Linux

COPY ./rtree /rtree
RUN R -e "setwd('/rtree');library(devtools);install_local('/rtree')"

RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager);BiocManager::install('rhdf5')"
RUN R -e "library(devtools);install_github('rformassspectrometry/ProtGenerics');install_github('rformassspectrometry/MsCoreUtils');install_github('rformassspectrometry/Spectra');install_github('rformassspectrometry/MsBackendMgf')"
RUN R -e "library(BiocManager);BiocManager::install('RcppArmadillo')"
RUN R -e "library(BiocManager);BiocManager::install('BiocParallel')"
RUN R -e "install.packages(c('XML','stringr','DBI','RSQLite','igraph','Matrix','jsonlite','optparse','tools'))"
RUN R -e "install.packages('data.table')"

RUN R -e "library(BiocManager);BiocManager::install('mzR')"
RUN R -e "library(BiocManager);BiocManager::install('MSnbase')"
RUN R -e "library(BiocManager);BiocManager::install('xcms')"
RUN R -e "library(BiocManager);BiocManager::install('cliqueMS')"
RUN R -e "library(BiocManager);BiocManager::install('CAMERA')"
RUN R -e "library(BiocManager);BiocManager::install('InterpretMSSpectrum')"

COPY onlineLCMSaligner /onlineLCMSaligner
#RUN apt-get install -y libgmp3-dev
RUN R -e "install.packages(c('ClusterR','ggplot2','gghighlight','lpSolve'))"
RUN R -e "setwd('/onlineLCMSaligner');library(devtools);install_local('/onlineLCMSaligner')"

COPY MZmineXMLManipulator /MZmineXMLManipulator
RUN R -e "setwd('/MZmineXMLManipulator');library(devtools);install_local('/MZmineXMLManipulator')"

COPY pylcmsprocessing /pylcmsprocessing

#This is the 2 workflow running script
COPY run_lcms_processing.sh /run_workflow.sh
COPY wrapper_docker.py /wrapper_docker.py
ENV LC_ALL C
ENV PATH="/:${PATH}"
ENTRYPOINT bash /run_workflow.sh
