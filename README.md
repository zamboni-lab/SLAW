# Zamboni Lab LCMSprocessing pipeline

## Current version : 0.9

The LCMS processing which has been developped has mainly been built for scalability and
efficiency, meaning that it supposed to be able to process thousands of samples and it should
extract them in a sensitive and reliable manner. However it was not design for user-friendliness.
If you want to manually inspect some small dataset, consider using compoundsDiscoverer or MZmine.

Therefore it is composed of different piece of software wrapped together, this to avoid installation, all the softwares are bundled together in a *Docker*. A docker can be seen as a mini virtual machine.

## Installing docker

Docker can be installed following this [tutorial](https://runnable.com/docker/install-docker-on-windows-10).
The whole installation process should take less than 15 minutes, mainly because you need to restart your computer.

After the installation is finished you have to add yourself to the docker-user group of ETH to run Docker. At the moment you have to do it yourself, to do so sign out of your session (Important : sign out) and log in the admin account of your computer (SYSBCPU) and find the local users group in the control panel :
![alt text](imgs/local_users.png)

Then on the windows click on the Groups folder in the middle panel, open the docker-users group. You  can then add yourself to the docker-users, using you ETH login (dalexis in my case). You can then log out of the SYSBCPU account
### Increasing processing power.
By default the docker virtual machine only take a small part of your computing power, however especially if you have a workstation. To do so right click on the docker icon at the right of your task bar :
 ![alt text](imgs/icon_docker.png)

And open the settings page :
 ![alt text](imgs/increasing_power.png)

The number of CPUs and the RAM will notably increase the speed of the peakpicking. One peakpicking experiment take 1.25Go of RAM approximately. After these modifications the docker engine will restart, which can take several minutes.

## Running the workflow

### Getting the last version of the workflow
Once docker is installed please run the Windows Powershell on windows or any Unix terminal. The workflow can then be installed directly form the "dockerhub" just by typing the following lines :
```
docker pull adelabriere/lcms_workflow_zamboni:latest
```
The required space is approximately 3 Go.

### Running the docker

_TLDR :_
```
docker run -it --cap-add=SYS_ADMIN --cap-add=DAC_READ_SEARCH --privileged -e INPUT=/sauer1/users/Alexis/examples_lcms_workflow/input -e OUTPUT=/sauer1/users/Alexis/examples_lcms_workflow/output -e USERNAME=dalexis adelabriere/lcms_workflow_zamboni:latest
```
with :
* __INPUT__ the input directory, which can be on Sauer1, here
* __OUTPUT__ the output directory which can also be on Sauer1, it should be an empty directory. It can potentially include a *parameters.txt* file
* __USERNAME__ your username to access sauer1 (the same than your windows session login)


###More detail

### What's inside
The docker is focused on the passage form raw-files to a as set of usable informations. The LC-MS workflow extract 2 kind informations, pekas and ms-ms spectra specific of MS-MS spectra.
* A set of peaktable containing all the peaks extracted in the MS1 format one for each .mzML files in a csv format each including the following fields :
  * _mz_ : The measured mass-to-charge ratio measured as the intensity weighted mz across the different scans.
  * _rt_ : The measured retention time at the apex of the peak
  * _height_ : The peak height in count or as an intensity
  * _intensity_ : The integrated area of the peak in the time dimension
  * _rt_min_ : The detected start of the peak in minute
  * _rt_max_ : The detected end of the peak in minute
  * _mz_min_ : The detected end of the peak in minute
  * _mz_max_ : The detected end of the peak in minute
  * _SN_ : The signal-to-noise ratio as measured by the ADAP algorithm
  * _peakwidth_ : The peakwidth of the chromatographic peak calculated as _rt_max - rt_min_
  * _right_on_left_assymetry_ : A measure of the peak assymetry calculated as  _(rt_max - rt)/(rt - rt_min)_
* A set of _.mgf_ containing all the MS-MS spectra in .mgf format. The .mgf format is standard format used by many tools to idetify ms-ms spectra (GNPS spectral search, SIRIUS, CFM-ID).


```
docker run -it --cap-add=SYS_ADMIN --cap-add=DAC_READ_SEARCH --privileged -e INPUT=/sauer1/users/Alexis/examples_lcms_workflow/input -e OUTPUT=/sauer1/users/Alexis/examples_lcms_workflow/output -e USERNAME=dalexis adelabriere/lcms_workflow_zamboni:latest
```


## What is inside this docker ?

The LCMS processing workflow incorporates three main steps :
- Peak picking using the MZmine workflow
- Peak alignment using an in-house aligner
- Peak feature annotations using the MScliques worklow.
- The MS-MS spectra are eventually output to be the majority of the daya.
