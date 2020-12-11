###Deterine the polarity of file

args <- commandArgs(trailingOnly = TRUE)

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
