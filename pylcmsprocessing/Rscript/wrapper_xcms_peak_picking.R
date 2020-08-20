suppressWarnings(suppressMessages(library(xcms, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))

args <- commandArgs(trailingOnly = TRUE)

# print(args)
# args <- c(
#   "U:/users/Alexis/examples_lcms_workflow/input/BEH30_2min_LipidMix_3.mzML",
#   "U:/users/Alexis/examples_lcms_workflow/BEH30_2min_LipidMix_3.csv",
#   "15",
#   "2",
#   "10",
#   "5",
#   "3",
#   "5000",
#   "5000"
# )

PATH_RAW <- args[1]
PATH_OUTPUT <- args[2]
PPM <- as.numeric(args[3])
PEAKWDITH <- as.numeric(c(args[4],args[5]))
SNT <- as.numeric(args[6])
PREFILTER <- as.numeric(c(args[7],args[8]))
NOISE <- as.numeric(args[9])


xraw <- xcmsRaw(PATH_RAW)
peaks <- xcms:::findPeaks.centWave(xraw, ppm = PPM, peakwidth = PEAKWDITH, 
  snthresh = SNT, prefilter = PREFILTER) 


dtable <- data.table(mz=peaks[,"mz"],rt=peaks[,"rt"]/60,height=peaks[,"maxo"],intensity=peaks[,"into"],
                     rt_min=peaks[,"rtmin"]/60,rt_max=peaks[,"rtmax"]/60,
                     mz_min=peaks[,"mzmin"],mz_max=peaks[,"mzmax"],SN=peaks[,"sn"],
                     peakwidth=(peaks[,"rtmax"]-peaks[,"rtmin"])/60)

fwrite(dtable,file = PATH_OUTPUT,col.names = TRUE)
