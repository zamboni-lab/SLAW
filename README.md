# SLAW

### Recent changes:
__18/11/2021__: Memory efficient spectral merging. The memory efficiency of ms-ms merging has been upscaled following bug report. It is still in dev as I am running some tests, but if you had some issues running SLAW with an OOM error after the datamatrix was generated please try with adelabriere/slaw:dev. 

### Introduction

SLAW is a scalable, containerized workflow for untargeted LC-MS processing. It was developed by Alexis Delabriere in the [Zamboni Lab](https://imsb.ethz.ch/research/zamboni.html) at ETH Zurich. An explanation of the advantages of SLAW and its motivations of development can be found in this [blog post](https://metabolomics.blog/2021/07/02/slaw/). In brief, the core advantages of SLAW are:
 * Complete processing including peak picking, sample alignment, pick picking, grouping of isotopologues and adducts, gap-filling by data recursion, extraction of consolidated MS2 spectra and isotopic data.
 * Scalability: SLAW can process thousands of samples efficiently
 * Wrapping of three main peak picking algorithms: Centwave, FeatureFinderMetabo, ADAP
 * Automated parameter optimization for picking, alignment, gap-filling

If you want to use SLAW, please cite the following paper:

> Delabriere A, Warmer P, Brennsteiner V and Zamboni N, *SLAW: A scalable and self-optimizing processing workflow for untargeted LC-MS*, Anal. Chem. 2021 (https://doi.org/10.1021/acs.analchem.1c02687)
 
### This repository contains the current stable version (1.0.0)**

The latest development version can be found on https://github.com/zamboni-lab/SLAW. 

Additional information is available on the paper and in the wiki section.

## Installation

The source code provided here is meant for developers. For users, setting up an environment with R, python, mzMine, etc. is a cumbersome process. Instead, **the recommended way to use SLAW is to pull the container from *DockerHub***. The SLAW container comes preconfigured with all components and - thanks to the self-optimizing algorithms - it can be used as a black box:

```
docker pull adelabriere/slaw:latest
```

If you need to run SLAW on a HPC infrastructure that does not provide the rights to operate Docker, an equivalent container is available from Singularity. As SingularityHub went read-only, the recommended way to get the singularity container is currently to pull it from DockerHub:

```
 singularity pull slaw.sif docker://adelabriere/slaw:latest
```

## Quick start

In principle, SLAW is executed by running the Docker container in terminal window (Linux) or in Powershell (Windows, don't use Powershell ISE). The minimum requirements are an input folder with **centroided** mzML files, and an output folder. Both are mounted with the -v option to either /input and /output: 
```
docker run --rm -v MZML_FOLDER:/input -v PATH_OUTPUT:/output adelabriere/slaw:latest
```
_Note: we recommend using a local folder as output. Network mounts can create problems._

_Note: centroided mzML can be obtained with ProteoWizard. Discard profile data and always prioritize the centroid data provided by vendor's software over recalculating centroids with pwiz. (Filter peakPicking vendor msLevel=1-2)_

_Note: more CPU cores will accelerate computation. If you are running parallel containers, the number can be set with -e NCORES=xy. However a high number of cores (>80) may encounter internal R limitations, there keep this number below 80._

_Note: a verbose execution can be activated with -e LOGGING=DEBUG_

Therefore, a more complete command-line example on a Windows machine is:
```
docker run --rm -v D:\mydata\input_folder_with_mzML:/input -v D:\mydata\output_folder:/output -e NCORES=16 -e LOGGING=DEBUG adelabriere/slaw:latest
```

If you specified the path correctly, you should see the following text:
```
2020-12-02|12:39:28|INFO: Total memory available: 7556 and 6 cores. The workflow will use 1257 Mb by core on 5 cores.
....
2020-12-02|12:39:31|INFO: Parameters file generated please check the parameters values/ranges and toggle optimization if needed.
```
If the output folder did not contain a _parameters.txt_ file (as in the previous example), SLAW crates a default version and stops. This allows to adjust settings. Upon editing, processing can then be resumed by running the same command as above. 

If the _parameters.txt_ file existed, SLAW reads all settings and proceeds with execution. 
```
2020-12-02|12:39:37|INFO: Total memory available: 7553 and 6 cores. The workflow will use 1257 Mb by core on 5 cores.
2020-12-02|12:39:37|INFO: Guessing polarity from file:DDA1.mzML
2020-12-02|12:39:38|INFO: Polarity detected: positive
2020-12-02|12:39:39|INFO: STEP: initialisation TOTAL_TIME:2.05s LAST_STEP:2.05s
...
2020-12-02|12:41:04|INFO: STEP: gap-filling TOTAL_TIME:86.86s LAST_STEP:15.17s
2020-12-02|12:41:30|INFO: Annotation finished
2020-12-02|12:41:30|INFO: STEP: annotation TOTAL_TIME:112.74s LAST_STEP:25.87s
```
A similar processing can be run using singularity like this:
```
singularity run -C -W workspace_folder -B PATH_OUTPUT:/output  -B MZML_FOLDER:/input slaw.sif
```
Some example data are given in the test_data folder of this repository. These data have been heavily filtered to allow quick testing (i.e. we removed low abundant centroids). An example of input folder is given in the _test_data/mzML_ folder, and an example of parameters file (which will be generated by SLAW if you run it on an empty folder) is given in _test_data/parameters.txt_ an example of the complete output of SLAW without optimization is given in _test_data/output_. 

Output files are generated in PATH_OUTPUT and include:
 * datamatrices: The complete table with rows corresponding to features (ions) and the columns corresponding to samples. Three types of matrices are generated.
 * fused_mgf: The consensus MS2 spectra consolidated in a single MGF (Mascot Generic Format) file. For each feature in the data matrix we report MS2 spectra for each collisional energy (if available) and a fused one.
 * OPENMS/CENTWAVE/ADAP: Stores the individual peak tables and ms-ms spectra for each sample as obtained from the peak pickers.
 * temp_optim: if optimization was done, this folder includes all key files and results created during parameter optimization

## Parameters
The _parameters.txt_ includes all settings used in processing, organized by topic. More details are provided in the file. In most cases, the _parameters.txt_ file lists a value and a min-max range. The range is only used when SLAW is asked to optimize parameters (see paper for details).

The most important parameter is **optimization > need_optimization > true/false**, and indicates whether SLAW will attemp to optimize parameters (i.e. those flagged with _essential_) using QC files. Depending on the size of the files and the peak picking algorithm, optimization can take several hours. To avoid long processing time, the optimization is switched off (False) by default. However, the default settings are unlikely to work for all LC-MS settings. Therefore, it is surely worth to execute it at least once by switchinig the parameter to True. The optimization will produced a _parameters.txt_ which is optimal for the given data set. If the chromatographic and MS conditions are stable, the optimized _parameters.txt_ file can simply be reused in subsequent runs with optmization toggled off (False). 

Peak picking: three algorithms are included: centWave (XCMS), ADAP (mzMine), and FeatureFinderMetabo (openMS). Depending on the type of instrument, scan rate, baseline, peak shapes,... they perform differently. In our experience, centWave is generally quite robust and our default. FeatureFinderMetabo in openMS is faster, and therefore of particular interest for large studies, large files, or for optimizing parameters related to chromtographic peak shape.
