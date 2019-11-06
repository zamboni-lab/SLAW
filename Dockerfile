FROM basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

COPY pylcmsprocessing /pylcmsprocessing

#Workflow running script
COPY run_lcms_processing.sh /run_workflow.sh

#We ocpy tht eMZmine software
COPY MZmine-2.52-Linux /MZmine-2.52-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

COPY onlineLCMSaligner /onlineLCMSaligner

RUN R -e "library(devtools);remove.packages('onlineLCMSaligner');install('/onlineLCMSaligner')"

RUN pip3 install --no-cache-dir wheel psutil


ENTRYPOINT bash run_workflow.sh
