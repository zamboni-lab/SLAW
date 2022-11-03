FROM zambonilab/slaw_basis
# Base image should always include all the dependencies

# THis make the DOckerfile build work on github CI
RUN apt-get update
RUN apt-get -y --no-install-recommends --fix-missing install openms

# Fix for singularity on HPC
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager);BiocManager::install('rhdf5')"
RUN R -e "library(devtools);install_github('rformassspectrometry/ProtGenerics');install_github('rformassspectrometry/MsCoreUtils');install_github('rformassspectrometry/Spectra');install_github('rformassspectrometry/MsBackendMgf')"

#Resinstalling data.table as it seems to become problematic after Rhdf5
RUN R -e "remove.packages('data.table');install.packages('data.table')"
RUN R -e "install.packages(c('XML','stringr','DBI','RSQLite','igraph','Matrix','jsonlite','optparse','tools'))"
RUN R -e "library(BiocManager);BiocManager::install('RcppArmadillo')"
RUN R -e "library(BiocManager);BiocManager::install('BiocParallel')"

# adding the dev version of netcdf
RUN apt-get install -y libnetcdf-dev
RUN R -e "library(BiocManager);BiocManager::install('mzR')"
RUN R -e "library(BiocManager);BiocManager::install('MSnbase')"
RUN R -e "library(BiocManager);BiocManager::install('xcms')"
RUN R -e "library(BiocManager);BiocManager::install('cliqueMS')"
RUN R -e "library(BiocManager);BiocManager::install('CAMERA')"
RUN R -e "library(BiocManager);BiocManager::install('InterpretMSSpectrum')"

##We install the 2 packages in the same folder
COPY onlineLCMSaligner /onlineLCMSaligner
RUN apt-get install -y libgmp3-dev
RUN R -e "install.packages(c('ClusterR',''ggplot2','gghighlight'))"
RUN R -e "setwd('/onlineLCMSaligner');library(devtools);install_local('/onlineLCMSaligner')"

COPY MZmineXMLManipulator /MZmineXMLManipulator
RUN R -e "setwd('/MZmineXMLManipulator');library(devtools);install_local('/MZmineXMLManipulator')"

#Dependency copy
COPY MZmine-2.52-Linux /MZmine-2.52-Linux
COPY pylcmsprocessing /pylcmsprocessing

#This is the 2 workflow running script
COPY run_lcms_processing.sh /run_workflow.sh
COPY wrapper_docker.py /wrapper_docker.py
ENV LC_ALL C
ENV PATH="/:${PATH}"
ENTRYPOINT bash /run_workflow.sh
