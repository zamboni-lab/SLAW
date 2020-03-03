FROM adelabriere/basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

RUN pip3 uninstall -y numpy
RUN pip3 install numpy==1.17


##The copying is always a file.
COPY onlineLCMSaligner /onlineLCMSaligner
RUN R -e "setwd('/onlineLCMSaligner');library(devtools);remove.packages('onlineLCMSaligner');install(pkg='/onlineLCMSaligner')"

COPY pylcmsprocessing /pylcmsprocessing
#Workflow running script
COPY run_lcms_processing.sh /run_workflow.sh
#We ocpy tht eMZmine software
COPY MZmine-2.52-Linux /MZmine-2.52-Linux
#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py


ENTRYPOINT bash run_workflow.sh
