#FROM adelabriere/basis_workflow
FROM zambonilab/slaw_newbase
# Base image should always include all the dependencies

# THis make the DOckerfile build work on github CI
RUN apt-get update
RUN apt-get -y --no-install-recommends --fix-missing install openms

# Fix for singularity on HPC
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

##We install the 3 packages in the same folder
#RUN R -e "install.packages('rtree')"
#COPY onlineLCMSaligner /onlineLCMSaligner
#RUN R -e "setwd('/onlineLCMSaligner');library(devtools);install_local('/onlineLCMSaligner')"
#COPY MZmineXMLManipulator /MZmineXMLManipulator
#RUN R -e "setwd('/MZmineXMLManipulator');library(devtools);install_local('/MZmineXMLManipulator')"
RUN R -e "library(BiocManager);BiocManager::install('rhdf5')"
RUN R -e "library(devtools);install_github('rformassspectrometry/ProtGenerics');install_github('rformassspectrometry/MsCoreUtils');install_github('rformassspectrometry/Spectra');install_github('rformassspectrometry/MsBackendMgf')"

#Resinstalling data.table as it seems to become problematic after Rhdf5
RUN R -e "remove.packages('data.table');install.packages('data.table')"

#Dependency copy
#COPY MZmine-2.52-Linux /MZmine-2.52-Linux
COPY pylcmsprocessing /pylcmsprocessing

#This is the 2 workflow running script
COPY run_lcms_processing.sh /run_workflow.sh
COPY wrapper_docker.py /wrapper_docker.py
ENV LC_ALL C
ENV PATH="/:${PATH}"
ENTRYPOINT bash /run_workflow.sh
