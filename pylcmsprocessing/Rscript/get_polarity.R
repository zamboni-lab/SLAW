###Deterine the polarity of file

args <- commandArgs(trailingOnly = TRUE)

# args <- c("U:/users/Alexis/examples_lcms_workflow/input/BEH30_2min_LipidMix_3.mzML","U:/users/Alexis/examples_lcms_workflow/pol.txt")

PATH_RAW <- args[1]
OUTPUT <- args[2]
suppressMessages(suppressWarnings(library(mzR,quietly = TRUE)))
xraw <- openMSfile(PATH_RAW)
pol <- unique(header(xraw)[,"polarity"])
if( pol==1 || pol=="1" || startsWith(tolower(pol),"pos")){
  pol <- "positive"
}else{
  pol <- "negative"
}

f <- file(OUTPUT,"w")
writeLines(pol,con=f)
close(f)
