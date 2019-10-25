FROM adelabriere/basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

#We copy the MZmine software
COPY pylcmsprocessing /pylcmsprocessing
COPY MZmineXMLManipulator /MZmineXMLManipulator

#Workflow running script
COPY run_lcms_processing.sh /run_workflow.sh

COPY MZmine-2.52-Linux /MZmine-2.52-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

RUN R -e "remove.packages('MZmineXMLManipulator');install.packages('/MZmineXMLManipulator',repos=NULL,type='source')"
#We install the mounting the dependencies
RUN pip3 install psutil

ENTRYPOINT bash run_workflow.sh
