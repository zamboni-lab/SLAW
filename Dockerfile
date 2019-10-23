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

COPY MZmine-2.51-Linux /MZmine-2.51-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

#We install the mounting the dependencies
RUN mkdir sauer1

ENTRYPOINT mount //nas22.ethz.ch/biol_imsb_sauer_1 sauer1 -nobrl -rw -o domain=d.ethz.ch,username=$USERNAME,vers=3.0 && bash
