# SLAW

## Note to users

This (public) repository  is updated and maintained according to users' requests and reported issues. SLAW is under continuous and active development in our lab. As this work includes a lot of customization to our workflows and environment, we don't publish novel features unless there is a specific request. Hence, please use the Issues tab to submit your ideas or problems.   

## Changelog

Can be found in the [NEWS](NEWS.md).

## Introduction

SLAW is a scalable, containerized workflow for untargeted LC-MS processing (DDA only).
It was developed by Alexis Delabriere in the [Zamboni Lab](https://imsb.ethz.ch/research/zamboni.html) at ETH Zurich.
An explanation of the advantages of SLAW and its motivations of development can be found in this [blog post](https://metabolomics.blog/2021/07/02/slaw/).
In brief, the core advantages of SLAW are:

* Complete processing including peak picking, sample alignment, pick picking, grouping of isotopologues and adducts, gap-filling by data recursion, extraction of consolidated MS2 spectra (only for DDA experiments! DIA-MS2 spectra will be skipped) and isotopic data
* Scalability: SLAW can process thousands of samples efficiently
* Wrapping of three main peak picking algorithms: Centwave, FeatureFinderMetabo, ADAP
* Automated parameter optimization for picking, alignment, gap-filling

If you want to use SLAW, please cite the following paper:

> Delabriere A, Warmer P, Brennsteiner V and Zamboni N, *SLAW: A scalable and self-optimizing processing workflow for untargeted LC-MS*, Anal. Chem. 2021 (https://doi.org/10.1021/acs.analchem.1c02687)

## Technical description

Please refer to the paper.
Additional information is given in the wiki section.

## Installation

The source code provided here is meant for developers.
For users, **the recommended way to use SLAW is to pull the container from *DockerHub***.
The SLAW container comes preconfigured with all components and - thanks to the self-optimizing algorithms - it can be used as a black box:

```
docker pull zambonilab/slaw:latest
```

## Input files

SLAW uses three (types of) inputs that have to be in a **unique folder**:

* **Raw MS data in mzML format**.
Files can include MS1 and DDA-MS2 scans.
DIA-MS is not supported.
**All data must be centroided and of unique polarity**.
Centroided mzML can be obtained with ProteoWizard.
Discard all profile data and always prioritize the centroid data provided by vendor's software, e.g. with the filter `peakPicking vendor msLevel=1-2`.
* A `parameters.txt` file including all preferences and settings that govern SLAW's behavior.
If a `parameters.txt` is missing, a default file is created upon executing the Docker.
In this case, execution of SLAW stops.
More details on parameters are provided below. 
* A `samples.csv` file that defines the file types.
In principle, this csv file consists of two columns (no headers).
The first column (from the left) reports the file names (with mzML ending).
Only the basename to the sample will be considered, e.g. if the name includes a full path like `C:/some/path/to/sample1.mzML` only the `sample1.mzML` part will be matched against the file names found in the same folder.
From top to bottom, rows should follow the order of injection.
The second column (from the left) must be one of the following four (case insensitive):
  * **QC**: QC samples used as reference during parameter optimization, to extract reference peaks for RT alignment, to detect isotopes/adducts, etc.
  Ideally, QC samples are pooled study samples that are scattered (intercalated) through the whole sequences.
  If no QC samples are defined, or the `samples.csv` is missing, SLAW will pick random sample files.
  * **samples**: Samples of all kind (study samples, calibrants, spike-ins, etc.).
  They should contain MS1 scans and optionally MS2 data. 
  * **MS2**: indicates files that include primarily MS2 spectra and should not be considered for MS1 peak picking and quantification.
  MS2 spectra will be mapped to MS1 features (identified in QC and samples) after alignment.
  MS2 is typically used to flag files used with targeted MS2, iterative MS2, or generally DDA-heavy files.
  * **blank**: Blanks will be discarded from the final peaktable.

A well-formatted example of `*.csv` file is (from test_data):

```
path,type
QC1.mzML,QC
S1.mzML,sample
S2.mzML,sample
DDA1.mzML,MS2
DDA2.mzML,MS2
QC2.mzML,QC
```

The `parameters.txt` file is mandatory.
If missing, SLAW will create a default file and stop.
Upon editing, you'll have to restart processing.

The `*.csv` file is optional.
If missing, SLAW will pick random samples and proceed with optimization.
  
## Quick start

SLAW is executed by running the Docker container in  a terminal window (on Linux) or in the Powershell (on Windows, don't use Powershell ISE).
All input files (all mzML files and, optionally `parameters.txt` and `samples.csv`) must be present in a single input folder.
In addition, **a local output folder must be created**.
Both folders must be mounted with the `-v` option to either `/input` and `/output` when running the container:

```
docker run --rm -v PATH_INPUT:/input -v PATH_OUTPUT:/output zambonilab/slaw:latest
```

In the above example, `PATH_INPUT` and `PATH_OUTPUT` should be replaced with the full path of the aforementioned input and output folder.
An example can be found at the end of this section.
Copying without updating will trigger an error because the paths are invalid. 

If you specified the paths correctly and they exist, you should see the following text:

```
2020-12-02|12:39:28|INFO: Total memory available: 7556 and 6 cores. The workflow will use 1257 Mb by core on 5 cores.
....
2020-12-02|12:39:31|INFO: Parameters file generated please check the parameters values/ranges and toggle optimization if needed.
```

**If the output folder did not contain a `parameters.txt` file**, SLAW creates a default version and stops.
This allows to adjust settings.
Upon editing, processing can then be resumed by running the same command as above.

**If the `parameters.txt` file existed**, SLAW reads all settings and proceeds with execution.

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

Example data are provided in the `test_data` folder of this repository.
These data have been heavily filtered to allow quick testing (i.e. we removed low abundant centroids).
An example of input folder is given in the `test_data/mzML` folder, and an example of parameters file (which will be generated by SLAW if you run it on an empty folder) is given in `test_data/parameters.txt`.
 <!-- an example of the complete output of SLAW without optimization is given in `test_data/output`. -->

During execution, SLAW generates the following **output files** in `PATH_OUTPUT`:

* *datamatrices*:
The complete table with rows corresponding to features (ions) and the columns corresponding to samples.
Three types of matrices are generated: before clustering and annotation of adducts/ions, after, and after upon collapsing adducts and ions.
* *fused_mgf*:
The consensus MS2 spectra consolidated in a single MGF (Mascot Generic Format) file.
For each feature in the data matrix we report MS2 spectra for each collisional energy (if available) and a fused one.
* *OPENMS/CENTWAVE/ADAP*: 
Stores the individual peak tables and ms-ms spectra for each sample as obtained from the peak pickers.
* *temp_optim*:
If parameter optimization was done, this folder includes all key files and results created during parameter optimization.

### Notes

* We recommend using a local folder as output.
Network mounts create problems.
* More CPU cores will accelerate computation.
If you are running parallel containers, the number can be set with `-e NCORES=n`.
However a high number of cores (n>80) may encounter internal R limitations, there keep this number below `80`.
* A verbose execution can be activated with `-e LOGGING=DEBUG`.
Therefore, a more complete command-line example on a Windows machine is:

```
docker run --rm -v D:\mydata\input_folder_with_mzML:/input -v D:\mydata\output_folder:/output -e NCORES=16 -e LOGGING=DEBUG zambonilab/slaw:latest
```

## Parameters

The `parameters.txt` includes all settings used in processing, organized by topic.
More details are provided in the file.
In most cases, the `parameters.txt` file lists a value and a min-max range.
The range is only used when SLAW is asked to optimize parameters (see paper for details).

The most important parameter is `optimization > need_optimization > true/false`, and indicates whether SLAW will attemp to optimize parameters (i.e. those flagged with `essential`) using QC files.
Depending on the size of the files and the peak picking algorithm, optimization can take several hours.
To avoid long processing time, the optimization is switched off (`False`) by default.
However, the default settings are unlikely to work for all LC-MS settings.
Therefore, it is surely worth to execute it at least once by switchinig the parameter to `True`.
The optimization will produced a `parameters.txt` which is optimal for the given data set.
If the chromatographic and MS conditions are stable, the optimized `parameters.txt` file can simply be reused in subsequent runs with optmization toggled off (`False`).

The second most important parameter is `noise_level_MS1`.
It allows to cut low abundance centroids.
On TOF data, these features include artifacts (e.g. peak ringing) and noise that bias and slow down optimization.
By setting `noise_level_MS1` to a number between `20-500` for TOFs and `1E4-1E5` for Orbitraps, we cut the grass.
The threshold is instrument and scan-rate dependent.
The best way to define the number is to inspect visually a MS1 scan in the raw data.
FYR, we typically use 100-150 for Agilent QTOFs (6550, 6546) and Sciex 7600 ZenoTOF.  

Peak picking: three algorithms are included: centWave (XCMS), ADAP (mzMine), and FeatureFinderMetabo (openMS).
Depending on the type of instrument, scan rate, baseline, peak shapes, ..., they perform differently.
In our experience, centWave is generally quite robust and our default.
FeatureFinderMetabo in openMS is faster, and therefore of particular interest for large studies, large files, or for optimizing parameters related to chromtographic peak shape.

Additional information can be found in the Wiki section.

## Missing speed?

* Don't use ADAP with parameter optimization.
It's slower than the others, with obvious consequences on run time.
The fastest is openMS.
Based on our experience with different types of TOF and Orbi spectra (some of which is reported in the paper), we recommend to use centWave as default.
OpenMS might be preferred for fast processing of large LC-MS studies with thousands of injections.
* Increase `noise_level_MS1` to `10-100`.
* Use a local folder as input.

## Use on HPC nodes

If you need to run SLAW on a HPC infrastructure that does not provide the rights to operate Docker, an equivalent container is available from Singularity.
As SingularityHub went read-only, the recommended way to get the singularity container is to pull it from DockerHub:

```
 singularity pull slaw.sif docker://zambonilab/slaw:latest
```

A similar processing can be run using singularity like this:

```
singularity run -C -W workspace_folder -B PATH_OUTPUT:/output  -B MZML_FOLDER:/input slaw.sif
```
