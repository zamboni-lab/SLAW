FROM adelabriere/basis_workflow
#Base image should always include all the dependencies

###INPUT FOLDER :
#"/rawfiles"

###OUTPUT FOLDER :
#"/output

#We copy the MZmine software
COPY pylcmsprocessing /pylcmsprocessing

COPY MZmine-2.51-Linux /MZmine-2.51-Linux

#The data needs to be run inside the docker.
COPY wrapper_docker.py /wrapper_docker.py

#We run the workflow
CMD python3 wrapper_docker.py
