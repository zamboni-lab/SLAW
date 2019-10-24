FROM adelabriere/basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output
RUN apt-get install -y cifs-utils

#We copy the MZmine software
COPY pylcmsprocessing /pylcmsprocessing
COPY fstab_sauer1 /etc/fstab
COPY run_lcms_processing.sh /run_workflow.sh

COPY MZmine-2.51-Linux /MZmine-2.51-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

#We install the mounting the dependencies
RUN mkdir sauer1

ENTRYPOINT bash run_workflow.sh
