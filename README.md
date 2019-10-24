# LCMSprocessing pipeline.

The LCMS processing which has been developped has mainly been built for scalability and
efficiency, meaning that it supposed to be able to process thousands of samples and it should
extract them in a sensitive and reliable manner. However it was not design for user-friendliness.
If you want to manually inspect some small dataset, consider looking for compoundsDiscoverer or MZmine.

Therefore it is composed of different piece of software wrapped together, this to avoid installation, all the softwares are bundled together in a *Docker*.

A *Docker* is very similar to virtual machine from the user perspectives. In this documentations, two parts will be distinguished.

## Installing docker

Docker can be installed following this [tutorial](https://runnable.com/docker/install-docker-on-windows-10).
The whole installation process should take less than 5 minutes. You then need 2 things, first to add yourself to the docker-users group and increase the memory available to docker.

You then have to add yourself to the docker-user group of ETH to run Docker. Here is a detailed step by step tutorial.

First go to the local users group in the control panel :
![alt text](imgs/local_users.png)

Then on the windows click on the Groups folder in the middle panel, open the docker-users group. You  can then add yourself to the user, using you eth login.

### Increasing processing power.
By default the docker virtual machine only take a small part of your computing power, however especially if you have a workstation, it can be


## Setting up docker
Docker by default is very limited in processing power, therefore to speed up the processing you should increase the possible memory and cores which are available to docker.


## Running the docker

Now that Docker is installed you are ready to go, you can process your data using the run command.


 To do you need to specify two repository, one which contains
your .mzML file (.mzML is the standard format used in the hyphenated MS world, you can convert your file in .mzML using the MSconvert Proteowizard softare), and one empty repository which will store the data. If the processing stop or fail for any reasons external to the software (updating code) don't erase this folder, calculations can restart.

The main step are the following :
- 1 Pulling the Docker (only the first time)
- 2a You run the docker and a file paramters.yaml is created in the output directory.
- 2b You already have a patameters.yaml file, you copy it into the output directory
- 3 If needed you tune or modify the parameters with the HIGH or ESSENTIAL tag. The other can be let as it is.
- 4 Run the docker fetch your results in the output directory.

### Step 1: Pulling the docker
Pulling the docker can be accomplished.


### Step 1 Running the docker a first time:
The docker needs two things run :
* An input directory containing the data to be processed, a bunch of .mzML files.
* An output directory, in this directory the parameters will be placed aswell as the

At the run time you also need to input your session password to connect on sauer1 using your session.

The LCMS workflow can take as output data placed on the sauer1 *NAS* or on a local drive (C:, etc...). The sauer1 drive is by default mounted inside the docker under the path */sauer1*. A demo version on a small subset of data is included inside my sauer1 folder, to run it simply.

```
docker run -it --cap-add=SYS_ADMIN --cap-add=DAC_READ_SEARCH --privileged -e INPUT=/sauer1/users/Alexis/examples_lcms_workflow/input -e OUTPUT=/sauer1/users/Alexis/examples_lcms_workflow/output -e USERNAME=dalexis lcms_workflow_zamboni
```

### Step 3 : Running the docker.

The docker needs to be run with the following command.

docker run lcms_workflow_zamboni:latest -v ./input:/rawfiles/mzML -v ./output:/output

The docker run eventually the following ocmmand inside the following workflow


## What is inside this docker ?

The LCMS processing workflow incorporates three main steps :
- Peak picking using the MZmine workflow
- Peak alignment using an in-house aligner
- Peak feature annotations using the MScliques worklow.
- The MS-MS spectra are eventually output to be the majority of the daya.
