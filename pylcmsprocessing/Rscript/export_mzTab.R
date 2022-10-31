suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(data.table,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(MsBackendMgf,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))

args <- commandArgs(trailingOnly = TRUE)

DEBUG <- FALSE

if (DEBUG) {
  args <- c("D:\\SW\\SLAW_test_data_out\\temp_processing_db.sqlite",
            "D:\\SW\\SLAW_test_data_out\\data_741d552fefa0759df99c04af0d7f6562.mzTab",
            "plus")
}
##
PATH_DB <- args[1]
OUTPUT_FILE_PATH <- args[2]
MZTAB_FORMAT <- args[3]
INPUT_SLAW_FOLDER <- dirname(args[1])

## decide export format
APPENDMGF <- FALSE
if (grepl('.mztabplus',tolower(OUTPUT_FILE_PATH))) APPENDMGF <- TRUE
if (grepl('plus',tolower(MZTAB_FORMAT))) APPENDMGF <- TRUE

if (APPENDMGF) OUTPUT_FILE_PATH <- str_replace(OUTPUT_FILE_PATH,'.mzTab','.mzTabPlus')

SEPARATOR <- "\t" 

SOFTWARE_NAME <- "SLAW"

dbb <- dbConnect(RSQLite::SQLite(), PATH_DB)
HASH <- dbGetQuery(dbb, "SELECT hash FROM peakpicking")[,1]
all_peaktables <- dbGetQuery(dbb, "SELECT path,output_ms,types,level FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
polarity <- dbGetQuery(dbb, "SELECT polarity FROM common")[,1]
dbDisconnect(dbb)

mztab_id <- paste(SOFTWARE_NAME,"_",HASH,".mzTab",sep="")

mtd <- c(paste("MTD","mzTab-version","2.2.0-M",sep=SEPARATOR),
         paste("MTD","mzTab-id",mztab_id,sep=SEPARATOR),
         "",
         paste("MTD","cv[1]-label","MS",sep=SEPARATOR),
         paste("MTD","cv[1]-full_name","PSI-MS controlled vocabulary",sep=SEPARATOR),
         paste("MTD","cv[1]-version","4.1.67",sep=SEPARATOR),
         paste("MTD","cv[1]-uri","https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",sep=SEPARATOR),
         "",
         paste("MTD","small_molecule-quantification_unit","[MS, MS:1001843, MS1 feature maximum intensity, ]",sep=SEPARATOR),
         paste("MTD","small_molecule_feature-quantification_unit","[MS, MS:1001843, MS1 feature maximum intensity, ]",sep=SEPARATOR),
         paste("MTD","small_molecule-identification_reliability","[MS, MS:1002955, hr-ms compound identification confidence level, ]",sep=SEPARATOR),
         paste("MTD","id_confidence_measure[1]","[MS, MS:1002888, small molecule confidence measure, ]",sep=SEPARATOR),
         "",
         paste("MTD","software[1]","[MS, MS:1002205, ProteoWizard msconvert, unknown]",sep=SEPARATOR),
         paste("MTD","software[2]","[MS, MS:1002878, SLAW, 1]",sep=SEPARATOR),
         paste("MTD","quantification_method","[MS, MS:1001834, LC-MS label-free quantitation analysis, ]",sep=SEPARATOR),
         "",
         paste("MTD","database[1]","[, , \"no database\", null]",sep=SEPARATOR),
         paste("MTD","database[1]-prefix","null",sep=SEPARATOR),
         paste("MTD","database[1]-version","Unknown",sep=SEPARATOR),
         paste("MTD","database[1]-uri","null",sep=SEPARATOR),
         "")

#
# MTD	database[1]	[MIRIAM, MIR:00100079, HMDB, ]
# MTD	database[1]-prefix	hmdb
# MTD	database[1]-version	4.0
# MTD	database[1]-uri	http://www.hmdb.ca/
#   MTD	database[2]	[MIRIAM, MIR:00000002, CHEBI, ]
# MTD	database[2]-prefix	chebi
# MTD	database[2]-version	Unknown
# MTD	database[2]-uri	https://www.ebi.ac.uk/chebi/
#   MTD	database[3]	[MIRIAM, MIR:00000013, KEGG Compound, ]
# MTD	database[3]-prefix	c
# MTD	database[3]-version	Unknown
# MTD	database[3]-uri	https://www.genome.jp/kegg/compound/


### samples
for(i in 1:length(all_peaktables$path)){
  mtd <- c(mtd,
           paste("MTD",paste0("sample[",i,"]"),paste0("sample_",formatC(i, width=4, flag="0")),sep=SEPARATOR),
           paste("MTD",paste0("sample[",i,"]-description"),paste(all_peaktables$types[i],all_peaktables$level[i]),sep=SEPARATOR))
}
mtd <- c(mtd,"")

### add msrun

# str_id_format <- "[MS, MS:1000768, Thermo nativeid format, ]"
str_id_format <- "[MS, MS:1000777, spectrum identifier nativeid format, ]"

str_polarity <- "[MS, MS:1000130, positive scan, ]"
if(polarity=="positive"){
  str_polarity <- "[MS, MS:1000130, positive scan, ]"
}else{
  str_polarity <- "[MS, MS:1000129, negative scan, ]"
}

####TODO put the good keyword for file
str_format <- "[MS, MS:1000584, mzML file, ]"
format_file <- str_split(all_peaktables$path,pattern = fixed("."))[[1]][2]
if(format_file=="mzML"){
  str_format <- "[MS, MS:1000584, mzML file, ]"
}

# partial_mtd <- vector(mode="list",length = )
for(i in 1:length(all_peaktables$path)){
  mtd <- c(mtd,
           paste("MTD",paste0("ms_run[",i,"]-location"),basename(all_peaktables$path[i]),sep=SEPARATOR),
           paste("MTD",paste0("ms_run[",i,"]-format"),str_format,sep=SEPARATOR),
           paste("MTD",paste0("ms_run[",i,"]-id_format"),str_id_format,sep=SEPARATOR),
           paste("MTD",paste0("ms_run[",i,"]-scan_polarity[1]"),str_polarity,sep=SEPARATOR))
}
##We add the mgf produced if there is one
# PATH_MGF <- file.path(INPUT_SLAW_FOLDER,"fused_mgf")
# PATH_MGF <- list.files(PATH_MGF,full.names = TRUE)
# if(length(PATH_MGF)>0){
#   mtd <- c(mtd,
#            "",
#            paste("MTD",paste0("ms_run[",i+1,"]-location"),basename(PATH_MGF),sep=SEPARATOR),
#            paste("MTD",paste0("ms_run[",i+1,"]-format"),"[MS, MS:1001062, Mascot MGF file, ]",sep=SEPARATOR),
#            paste("MTD",paste0("ms_run[",i+1,"]-id_format"),"[MS, MS:1000774, multiple peak list nativeid format, ]",sep=SEPARATOR),
#            paste("MTD",paste0("ms_run[",i+1,"]-scan_polarity[1]"),str_polarity,sep=SEPARATOR)
#   )
#
# }
mtd <- c(mtd,"")

# assays

for(i in 1:length(all_peaktables$path)){
  mtd <- c(mtd,
           paste("MTD",paste0("assay[",i,"]"),paste0("seq_inj_",formatC(i, width=4, flag="0")),sep=SEPARATOR),
           paste("MTD",paste0("assay[",i,"]-sample_ref"),paste0("sample[",i,"]"),sep=SEPARATOR),
           paste("MTD",paste0("assay[",i,"]-ms_run_ref"),paste0("ms_run[",i,"]"),sep=SEPARATOR))
}
mtd <- c(mtd,"",
         paste("MTD","study_variable[1]","undefined",sep=SEPARATOR),
         paste("MTD","study_variable[1]_refs","undefined",sep=SEPARATOR),
         "")

com <- c(paste("COM",paste("File generated from folder:",INPUT_SLAW_FOLDER, "(",date(),")"),sep=SEPARATOR),"")

##All table are derived from the full table
# DIR_DM <- file.path(INPUT_SLAW_FOLDER,"datamatrices")
PATH_FULL <- list.files(INPUT_SLAW_FOLDER,"data_full.*.csv",full.names = TRUE)
if (length(PATH_FULL)>0) {
  ##SMF table Construction
  dm <- fread(PATH_FULL,sep="\t")
  if (dim(dm)[2]<2) dm <- fread(PATH_FULL,sep=";")
  if (dim(dm)[2]<2) dm <- fread(PATH_FULL,sep=",")
  ocnames <- colnames(dm)
  quant_prefix <- paste(str_split(ocnames[length(ocnames)],fixed("_"))[[1]][1],"_",sep="")
  quant_cols <- which(startsWith(ocnames,quant_prefix))

} else {
  ## use data_filled instead
  PATH_FULL <- list.files(file.path(INPUT_SLAW_FOLDER,"temp"),"data_filled*.csv",full.names = TRUE)
  ##SMF table Construction
  dm <- fread(PATH_FULL,sep="\t")
  if (dim(dm)[2]<2) dm <- fread(PATH_FULL,sep=";")
  if (dim(dm)[2]<2) dm <- fread(PATH_FULL,sep=",")
  ocnames <- colnames(dm)
  quant_prefix <- paste(str_split(ocnames[length(ocnames)],fixed("_"))[[1]][1],"_",sep="")
  quant_cols <- which(startsWith(ocnames,quant_prefix))

  # add missing columns
  dm$annotation <- rep("null",nrow(dm))
  dm$clique <- rep("null",nrow(dm))
  dm$group <- 1:nrow(dm)
}

SMF_ID <- 1:nrow(dm)
SME_ID_REFS <- rep("null",nrow(dm))
SME_ID_REF_ambiguity_code  <- rep("null",nrow(dm))

addduct_found <- table(dm$group)

count_clique <- dm$clique

treated_annot <- str_match(dm$annotation,"(\\[[0-9]?M[\\+\\-][a-zA-Z0-9\\+\\-]*\\])([0-9]?[\\+\\-]) ?\\+?([0-9])?")

smf_adduct_ion <- apply(treated_annot[,2:3],1, paste,  collapse= "")

####We put the sol annotation to null
non_identified <- seq_along(smf_adduct_ion)
num_evidence <- tapply(smf_adduct_ion,INDEX=as.factor(dm$group),FUN=function(x){length(unique(x))})
evidence_numbers <- num_evidence[match(as.character(dm$group),names(num_evidence))]
# smf_adduct_ion[evidence_numbers==1] <- "null"
smf_isotopomer <- treated_annot[,4]
smf_isotopomer[!(is.na(smf_isotopomer))] <- sapply(smf_isotopomer[!(is.na(smf_isotopomer))],
                                                   function(x){paste0('[MS, MS:1002956, isotopic ion MS peak, M+',x,' peak ]')})
# smf_isotopomer[(is.na(smf_isotopomer)) & (evidence_numbers>1)] <- 0
smf_isotopomer[(is.na(smf_isotopomer))] <- "null"
# smf_isotopomer[(smf_isotopomer==1)-1] <- 0 # risky!! but should do
opt_smf_isotopomer <- as.numeric(treated_annot[,4])
opt_smf_isotopomer[(is.na(opt_smf_isotopomer)) & (evidence_numbers>1)] <- 0
# opt_smf_isotopomer[(is.na(opt_smf_isotopomer)) & (evidence_numbers<2)] <- "null"
opt_smf_isotopomer[which(opt_smf_isotopomer==1)-1] <- 0 # risky!! but should do
opt_smf_isotopomer[is.na(opt_smf_isotopomer)] <- "null"
smf_exp_mass_to_charge <- dm$mz

opt_complete_ion <- dm$annotation

parseCharge <- function(x){
  if(x=="+") return(1)
  if(x=="-") return(-1)
  num <- as.numeric(str_sub(x,end = -2)[[1]])
  pol <- str_sub(x,start = 2)[[1]]
  if(pol=="-") return(-num)
  return(num)
}
smf_charge <- rep("null",nrow(dm))
known_charge <- !is.na(treated_annot[,3])
smf_charge[known_charge] <- sapply(treated_annot[known_charge,3],parseCharge)

smf_retention_time_in_seconds <- sprintf("%0.1f",dm$rt*60)
smf_retention_time_in_seconds_start <- sprintf("%0.1f",dm$min_rt_cor*60)
smf_retention_time_in_seconds_end <- sprintf("%0.1f",dm$max_rt_cor*60)

smf_abundance_assay <- ceiling(dm[,..quant_cols])
colnames(smf_abundance_assay) <- paste("abundance_assay[",1:length(quant_cols),"]",sep="")

smf_abundance_assay$'abundance_study_variable[1]' <- "null"
smf_abundance_assay$'abundance_variation_study_variable[1]' <- "null"


smf_opt_global_raw_isotopic_dist <- rep("null",nrow(dm))
smf_opt_global_raw_isotopic_name <- rep("null",nrow(dm))
if("isotopic_pattern_abs" %in% colnames(dm)){
  smf_opt_global_raw_isotopic_dist <- dm$isotopic_pattern_abs
  smf_opt_global_raw_isotopic_dist[smf_opt_global_raw_isotopic_dist==""] <- "null"
}
if("isotopic_pattern_annot" %in% colnames(dm)){
  smf_opt_global_raw_isotopic_name <- dm$isotopic_pattern_annot
  smf_opt_global_raw_isotopic_name[smf_opt_global_raw_isotopic_name==""] <- "null"
}
smf_opt_msms_index <- rep("null",nrow(dm))
smf_id_map <- data.frame()
if("mgf_ms2_id" %in% colnames(dm)){
  for (i in 1:nrow(dm)) {
    s <- dm$mgf_ms2_id[i]
    v <- suppressWarnings(na.omit(as.numeric(str_extract_all(s,"\\(?[0-9,.]+\\)?")[[1]])))
    if (length(v)>0) {
      smf_opt_msms_index[i] <- paste(v,collapse = '|')
      smf_id_map <- rbind(smf_id_map,data.frame(smf_id = replicate(length(v),SMF_ID[i]),mgf_id = v))
    }
  }
}

smf_opt_msms_index[is.na(smf_opt_msms_index)] <- "null"

smf_table <- data.table(SFH=rep("SMF",nrow(dm)),SMF_ID=SMF_ID, SME_ID_REFS = SME_ID_REFS,
                        SME_ID_REF_ambiguity_code= SME_ID_REF_ambiguity_code,
                        adduct_ion=smf_adduct_ion,
                        opt_global_ion = opt_complete_ion,
                        isotopomer=smf_isotopomer,
                        opt_global_isotopomer = opt_smf_isotopomer,
                        exp_mass_to_charge=smf_exp_mass_to_charge,charge=smf_charge,retention_time_in_seconds=smf_retention_time_in_seconds,
                        retention_time_in_seconds_start=smf_retention_time_in_seconds_start,
                        retention_time_in_seconds_end=smf_retention_time_in_seconds_end)

smf_table[["opt_global_raw_isotopic_dist"]] <- smf_opt_global_raw_isotopic_dist
smf_table[["opt_global_raw_isotopic_name"]] <- smf_opt_global_raw_isotopic_name
smf_table[["opt_global_mgf_index_ms2"]] <- smf_opt_msms_index

for(cc in colnames(smf_abundance_assay)){
  smf_table[[cc]] <- smf_abundance_assay[[cc]]
}

# SMH	SML_ID SMF_ID_REFS	database_identifier	chemical_formula	smiles	inchi	chemical_name	uri	theoretical_neutral_mass	adduct_ions	reliability	best_id_confidence_measure	best_id_confidence_value	abundance_assay[1]	abundance_assay[2]	abundance_assay[3]	abundance_assay[4]	abundance_assay[5]	abundance_assay[6]	abundance_study_variable[1]	abundance_variation_study_variable[1]	abundance_study_variable[2]	abundance_variation_study_variable[2]	opt_global_Progenesis_identifier
# SML	469	6 | 937	CHEBI:16737	C4H7N3O	null	null	Creatinine	null	113.0589	[M+H]+ | [M+Na]+	2	[MS,MS:1002889,Progenesis MetaScope score,]	56.4424	59809754.62	63773291.22	61630638.37	59131129.43	63874747.95	66320708.84	185213684.2	3.2135	189326586.2	5.7923	6.90_113.0582n
#

idx_abundance <- which(str_detect(colnames(smf_table),"abundance_assay"))
parse_subtable <- function(x,idx_abundance){
  sep_val <- "|"
  sml_smf_id_refs <- paste(x$SMF_ID,collapse = sep_val)
  sml_db_id <- "null"
  sml_chemical_formula <- "null"
  smiles <- "null"
  inchi<- "null"
  chemical_name<- "null"
  uri<- "null"
  theoretical_neutral_mass<- "null"
  # opt_adduct_ions <- paste(x$annotation,collapse = sep_val)
  adduct_ions <- paste(x$adduct_ion,collapse = sep_val)
  reliability <- "4"
  best_id_confidence_measure <- "null"
  abundances <- x[,..idx_abundance]
  cnames <- colnames(abundances)
  quant <- apply(abundances,2,sum)
  unique_ms2 <- unique(x$opt_global_mgf_index_ms2)
  ms2_str <- NULL
  if(unique_ms2[1]=="null"){
    ms2_str <- "null"
  }else{
    ms2_str <- paste(unique_ms2[unique_ms2!="null"],collapse = sep_val)
  }

  df <- list(SMF_ID_REFS=sml_smf_id_refs, database_identifier= sml_db_id,chemical_formula=sml_chemical_formula,smiles=smiles,inchi=inchi,
             chemical_name=chemical_name,uri=uri,theoretical_neutral_mass=theoretical_neutral_mass,adduct_ions=adduct_ions,
             reliability=reliability,best_id_confidence_measure=best_id_confidence_measure,best_id_confidence_value=best_id_confidence_measure,
             opt_global_mgf_index_ms2=ms2_str)
  df <- c(df,quant)
  data.frame(df,stringsAsFactors = FALSE,check.names = FALSE)
}

cnames <- colnames(smf_table)
sml_table <- by(cbind(smf_table,dm),INDICES = dm$group,FUN = parse_subtable,idx_abundance=idx_abundance)
sml_table <- do.call(rbind,sml_table)


sml_table$'abundance_study_variable[1]' <- "null"
sml_table$'abundance_variation_study_variable[1]' <- "null"

cnames <- colnames(sml_table)
sml_table <- cbind(rep("SML",nrow(sml_table)),1:nrow(sml_table),sml_table)
colnames(sml_table) <- c("SMH","SML_ID",cnames)


f <- file(OUTPUT_FILE_PATH,"w")
writeLines(com,f)
writeLines(mtd,f)
close(f)

fwrite(sml_table,OUTPUT_FILE_PATH,sep="\t",append = TRUE,col.names = TRUE)
f <- file(OUTPUT_FILE_PATH,"a")
writeLines(c(""),f)
close(f)
fwrite(smf_table,OUTPUT_FILE_PATH,sep="\t",append = TRUE,col.names = TRUE)

if (APPENDMGF==T) {
  MGF_DATA <- paste0(dirname(OUTPUT_FILE_PATH),"/spectra_",HASH,".mgf")
  if (!file.exists(MGF_DATA)) {return(NULL)}
  suppressWarnings(mgf <- readMgf(MGF_DATA))
  mgf_dat <- mgf@listData
  n <- mgf@nrows
  if ("acquisitionNum.1" %in% names(mgf_dat)) {mgf_dat$acquisitionNum <- mgf_dat$'acquisitionNum.1'}
  if (is.na(mgf_dat$acquisitionNum[1])) mgf_dat$acquisitionNum <- 1:n
  if ("PRECURSOR_INTENSITY" %in% names(mgf_dat)) {mgf_dat$precursorIntensity <- as.numeric(mgf_dat$'PRECURSOR_INTENSITY')}
  if ("ENERGY" %in% names(mgf_dat)) {mgf_dat$collisionEnergy <- as.numeric(mgf_dat$'ENERGY')}
  if ("MSLEVEL" %in% names(mgf_dat)) {mgf_dat$msLevel <- as.numeric(mgf_dat$'MSLEVEL')}
  df <- data.frame(COM = replicate(n,'COM'),
                   MGH = replicate(n,'MGF'),
                   mgf_id = as.numeric(mgf_dat$acquisitionNum),
                   prec_id = as.numeric(mgf_dat$FEAT_ID),
                   prec_rt = mgf_dat$rtime,
                   prec_mz = mgf_dat$precursorMz,
                   prec_int = mgf_dat$precursorIntensity,
                   energy = as.numeric(mgf_dat$collisionEnergy),
                   level = mgf_dat$msLevel,
                   title = mgf_dat$TITLE,
                   spec_mz = I(replicate(n,list())),
                   spec_int = I(replicate(n,list())),
                   spec_tic = replicate(n,NA),
                   spec_len = replicate(n,NA)
  )
  mz_end <- c(0,mgf_dat$mz@partitioning@end)
  int_end <- c(0,mgf_dat$int@partitioning@end)

  df <- merge(df,smf_id_map,by=c('mgf_id'),all.y = F)
  df$prec_id <- df$smf_id

  for (i in 1:nrow(df)){
    df[[i,'spec_mz']] <- mgf_dat$mz@unlistData[(mz_end[i]+1):mz_end[i+1]]
    df[[i,'spec_int']] <- mgf_dat$intensity@unlistData[(int_end[i]+1):int_end[i+1]]
    df[i,'spec_tic'] <- sum(df[[i,'spec_int']])
    df[i,'spec_len'] <- mz_end[i+1]-mz_end[i]
  }

  df <- df[order(df$mgf_id),c('COM','MGH','mgf_id','prec_id','prec_rt','prec_mz','prec_int','energy','level','title','spec_tic','spec_len','spec_mz','spec_int')]



  #df$MGH <- 'COM MGF'
  df[df==''] <- 'null'
  df[is.na(df)] <- 'null'
  f <- file(OUTPUT_FILE_PATH,"a")
  writeLines(c("","COM\tThe below part is NOT officially supported by the reference implementation of mzTag.",
               "COM\tIt's an appendix that includes all associated MGF data",""),f)
  close(f)
  fwrite(df,OUTPUT_FILE_PATH,sep="\t",append = TRUE,col.names = TRUE)

}
