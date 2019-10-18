# LCMSprocessing pipeline.

The LCMS processing which has been devlopped has mainly been built for scalability and
efficiency, meaning that it supposed to be able to process thousands of samples and it should
extract them in a sensitive and reliable manner. However it wa snot design for user-friendliness.
If you want to manually inspect some small dataset, I compounds discoverer aswelll as

Therefore it is composed of different piece of software wrapped together in different languages for
maximum speed and efficiency bundled in a *Docker*

The LCMS processing pipeline wrapped within a docker. A Docker is conceptually similar
to a virtual machine. All the dependencies aswell as the latest version of the pipeline have
bene installed on it, meaning that you can run the docker and that's it. We can discuss about the interactivity later.

Let's start with the basics :

## Installing docker
Docker can be installed following this [tutorial](https://runnable.com/docker/install-docker-on-windows-10).
The whole installation process should take less than 5 minutes.

You then have to add yourself to the docker-user group of ETH to run Docker. It is interesting because you can do it your self except if you
Here is a detailed step by step tutorial to do so :

TUTORIAL ADDING USERS

## Running the docker

Now that Docker is installed you are ready to go, you have to launch the docker
in order to process LC-MS pipeline. To do you need to specify two repository, one which contains
your .mzML file (.mzML is the standard format used in the hyphenated MS world, you can convert your file in .mzML using the MSconvert Proteowizard softare), and one empty repository which will store the data. If the processing stop or fail for any reasons external to the software (updating code) don t erase this repository, calculation can restart.

The main step are the following :
- 1a You run the docker and a file paramters.yaml is created in the output directory.
- 1b You already have a patameters.yaml file, you copy it into the output directory
- 2 If needed you tune or modify the parameters with the HIGH or ESSENTIAL tag. The other can be let as it is.
- 3 Run the docker fetch your results in the output directory.

### Step 3 : Running the docker.

The docker needs to be run with the following command.

docker run lcms_workflow_zamboni:latest -v ./input:/rawfiles/mzML -v ./output:/output

The docker run eventually the following ocmmand inside the following workflow


##What is inside this docker ?

The LCMS processing workflow incorporates three main steps :
- Peak picking using the MZmine workflow
- Peak alignment using an in-house aligner
- Peak feature annotations using the MScliques worklow.
- The MS-MS spectra are eventually output to be the majority of the daya.
