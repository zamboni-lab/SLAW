FROM adelabriere/basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

#We copy the MZmine software
RUN mkdir /temp_inst
COPY pylcmsprocessing /temp_inst/pylcmsprocessing
COPY MZmineXMLManipulator /temp_inst/MZmineXMLManipulator
COPY rtree /temp_inst/rtree
COPY onlineLCMSaligner /temp_inst/onlineLCMSaligner

RUN apt-get -y --no-install-recommends install libgmp3-dev
RUN R -e "library(devtools);install('/temp_inst/rtree');install.packages('lpSolve');install.packages('ClusterR');install.packages('viridis');remove.packages('MZmineXMLManipulator');install.packages('/temp_inst/MZmineXMLManipulator',repos=NULL,type='source');install('/temp_inst/onlineLCMSaligner')"
#Workflow running script
COPY run_lcms_processing.sh /run_workflow.sh

COPY MZmine-2.52-Linux /MZmine-2.52-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

#We install the mounting the dependencies
RUN pip3 install psutil

ENTRYPOINT bash run_workflow.sh
