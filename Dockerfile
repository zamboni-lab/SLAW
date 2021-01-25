FROM basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

#Moved in basis
#RUN pip3 uninstall -y numpy
#RUN pip3 install numpy==1.17
#RUN pip3 install sklearn statsmodels
RUN apt-get -y --no-install-recommends install openms

# Fix for singularity on HPC
RUN strip --remove-section=.note.ABI-tag /usr/lib/x86_64-linux-gnu/libQt5Core.so.5

##The copying is always a file.
COPY onlineLCMSaligner /onlineLCMSaligner
RUN R -e "setwd('/onlineLCMSaligner');library(devtools);remove.packages('onlineLCMSaligner');install(pkg='/onlineLCMSaligner')"
COPY MZmineXMLManipulator /MZmineXMLManipulator
RUN R -e "setwd('/MZmineXMLManipulator');library(devtools);remove.packages('MZmineXMLManipulator');install(pkg='/MZmineXMLManipulator')"

COPY pylcmsprocessing /pylcmsprocessing
#Workflow running script
COPY run_lcms_processing.sh /run_workflow.sh
#We ocpy tht eMZmine software
COPY MZmine-2.52-Linux /MZmine-2.52-Linux
#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py
ENTRYPOINT bash run_workflow.sh
