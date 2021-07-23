library(DBI)
library(RSQLite)
library(mzR)
library(stringr)
library(data.table)

###COnverting an output experiment folder into a mzTab format

args <- commandArgs(trailingOnly=TRUE)


###We create the header
EX_OUTPUT <- args[1]
PATH_OUTPUT <- args[2]
SEPARATOR <- "\t"

SOFTWARE_NAME <- "SLAW"
PATH_DB <- list.files(EX_OUTPUT,pattern=".sqlite",full.names = TRUE)

###If the temporary DB is also present
PATH_DB <- PATH_DB[length(PATH_DB)]

dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
id <- dbGetQuery(dbb, "SELECT hash FROM peakpicking")[,1]
###Getting the valid samples
all_peaktables <- dbGetQuery(dbb, "SELECT path,output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
polarity <- dbGetQuery(dbb, "SELECT polarity FROM common")[,1]
dbDisconnect(dbb)

mztab_id <- paste(SOFTWARE_NAME,"_",id,".mzTab",sep="")

####TODO put the good keyword for file
# f1 <- all_peaktables[1,1]
# f1 <- str_replace(f1,"/sauer1","U:")
# f1 <- openMSfile(f1)
str_id_format <- "[MS, MS:1000768, Thermo nativeID format, ]"


str_polarity <- "[MS, MS:1000130, positive scan, ]"
if(polarity=="positive"){
  str_polarity <- "[MS, MS:1000130, positive scan, ]"
}else{
  str_polarity <- "[MS, MS:1000129, negative scan, ]"
}

####TODO put the good keyword for file
format_file <- str_split(all_peaktables$path,pattern = fixed("."))[[1]][2]
str_format <- "[MS, MS:1000584, mzML file, ]"
if(format_file=="mzML"){
  str_format <- "[MS, MS:1000584, mzML file, ]"
}

partial_mtd <- vector(mode="list",length = nrow(all_peaktables))
for(i in seq_along(partial_mtd)){
  partial_mtd[[i]] <- c(
    paste(paste0("ms_run[",i,"]-location"),all_peaktables$path[i],sep=SEPARATOR),
    paste(paste0("ms_run[",i,"]-scan_polarity[1]"),str_polarity,sep=SEPARATOR),
    paste(paste0("ms_run[",i,"]-format"),str_format,sep=SEPARATOR),
    paste(paste0("ms_run[",i,"]-id_format"),str_id_format,sep=SEPARATOR)
  )
}

##We add the mgf produced if there is one
PATH_MGF <- file.path(EX_OUTPUT,"fused_mgf")
PATH_MGF <- list.files(PATH_MGF,full.names = TRUE)
if(length(PATH_MGF)>0){
  ms2_idx <- length(partial_mtd)+1
  partial_mtd[[ms2_idx]] <- c(
    paste(paste0("ms_run[",ms2_idx,"]-location"),PATH_MGF,sep=SEPARATOR),
    paste(paste0("ms_run[",ms2_idx,"]-scan_polarity[1]"),str_polarity,sep=SEPARATOR),
    paste(paste0("ms_run[",ms2_idx,"]-format"),"[MS, MS:1001062, Mascot MGF file, ]",sep=SEPARATOR),
    paste(paste0("ms_run[",ms2_idx,"]-id_format"),"[MS, MS:1000774, multiple peak list nativeID format, ]",sep=SEPARATOR)
  )

}


partial_mtd <- unlist(partial_mtd)

partial_mtd <- c(paste("mzTab-version","2.2.0-M",sep=SEPARATOR),
                 paste("mzTab-ID",mztab_id,sep=SEPARATOR),unlist(partial_mtd))

partial_mtd <- c(paste("MTD",partial_mtd,sep="\t"),"")


com_lines <- c(paste("File generated from folder:",EX_OUTPUT, "(",date(),")"))
com_lines <- paste("COM",com_lines,sep=SEPARATOR)

f <- file(PATH_OUTPUT,"w")
writeLines(com_lines,f)
writeLines(partial_mtd,f)
close(f)

##We get the annotated_peaktable
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
id <- dbGetQuery(dbb, "SELECT hash FROM peakpicking")[,1]
###Getting the valid samples
all_peaktables <- dbGetQuery(dbb, "SELECT path,output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
polarity <- dbGetQuery(dbb, "SELECT polarity FROM common")[,1]
dbDisconnect(dbb)


##All table are derived from the full table
DIR_DM <- file.path(EX_OUTPUT,"datamatrices")
PATH_FULL <- list.files(DIR_DM,"*full.csv",full.names = TRUE)


##SMF table Construction
dm <- fread(PATH_FULL,sep=";")
ocnames <- colnames(dm)
quant_prefix <- paste(str_split(ocnames[length(ocnames)],fixed("_"))[[1]][1],"_",sep="")
quant_cols <- which(startsWith(ocnames,quant_prefix))


SMF_ID <- 1:nrow(dm)

addduct_found <- table(dm$group)

count_clique <- dm$clique

### adduct_ion	isotopomer		charge are dervied form anntation

treated_annot <- str_match(dm$annotation,"(\\[[0-9]?M[\\+\\-][a-zA-Z0-9\\+\\-]*\\])([0-9]?[\\+\\-]) ?\\+?([0-9])?")

smf_adduct_ion <- apply(treated_annot[,2:3],1, paste,  collapse= "")

####We put the sol annotation to null
non_identified <- seq_along(smf_adduct_ion)
num_evidence <- tapply(smf_adduct_ion,INDEX=as.factor(dm$group),FUN=function(x){length(unique(x))})
evidence_numbers <- num_evidence[match(as.character(dm$group),names(num_evidence))]
smf_adduct_ion[evidence_numbers==1] <- "null"



smf_isotopomer <- treated_annot[,4]
smf_isotopomer[is.na(smf_isotopomer)] <- "null"
smf_exp_mass_to_charge <- dm$mz



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

smf_opt_global_raw_isotopic_dist <- rep("null",nrow(dm))
smf_opt_global_raw_isotopic_name <- rep("null",nrow(dm))
if("raw_isotopic_pattern" %in% colnames(dm)){
  smf_opt_global_raw_isotopic_dist <- dm$raw_isotopic_pattern
  smf_opt_global_raw_isotopic_dist[smf_opt_global_raw_isotopic_dist==""] <- "null"
}
if("raw_isotopic_pattern_annot" %in% colnames(dm)){
  smf_opt_global_raw_isotopic_name <- dm$isotopic_pattern_annot
  smf_opt_global_raw_isotopic_name[smf_opt_global_raw_isotopic_name==""] <- "null"
}
smf_opt_msms_index <- rep("null",nrow(dm))
if("ms2_id" %in% colnames(dm)){
  smf_opt_msms_index <- dm$ms2_id
  smf_opt_msms_index[is.na(smf_opt_msms_index)] <- "null"
}

smf_table <- data.table(SFH=rep("SMF",nrow(dm)),smf_id=SMF_ID,adduct_ion=smf_adduct_ion,isotopomer=smf_isotopomer,
                        exp_mass_to_charge=smf_exp_mass_to_charge,charge=smf_charge,retention_time_in_seconds=smf_retention_time_in_seconds,
                        retention_time_in_seconds_start=smf_retention_time_in_seconds_start,
                        retention_time_in_seconds_end=smf_retention_time_in_seconds_end)

for(cc in colnames(smf_abundance_assay)){
  smf_table[[cc]] <- smf_abundance_assay[[cc]]
}
smf_table[["opt_global_raw_isotopic_dist"]] <- smf_opt_global_raw_isotopic_dist
smf_table[["opt_global_raw_isotopic_name"]] <- smf_opt_global_raw_isotopic_name
smf_table[["opt_global_mgf_index_ms2"]] <- smf_opt_msms_index

# SMH	SML_ID	SMF_ID_REFS	database_identifier	chemical_formula	smiles	inchi	chemical_name	uri	theoretical_neutral_mass	adduct_ions	reliability	best_id_confidence_measure	best_id_confidence_value	abundance_assay[1]	abundance_assay[2]	abundance_assay[3]	abundance_assay[4]	abundance_assay[5]	abundance_assay[6]	abundance_study_variable[1]	abundance_variation_study_variable[1]	abundance_study_variable[2]	abundance_variation_study_variable[2]	opt_global_Progenesis_identifier
# SML	469	6 | 937	CHEBI:16737	C4H7N3O	null	null	Creatinine	null	113.0589	[M+H]+ | [M+Na]+	2	[MS,MS:1002889,Progenesis MetaScope score,]	56.4424	59809754.62	63773291.22	61630638.37	59131129.43	63874747.95	66320708.84	185213684.2	3.2135	189326586.2	5.7923	6.90_113.0582n
#

idx_abundance <- which(str_detect(colnames(smf_table),"abundance_assay"))
parse_subtable <- function(x,idx_abundance){
  sep_val <- "|"
  sml_smf_id_refs <- paste(x$smf_id,collapse = sep_val)
  sml_db_id <- "null"
  sml_chemical_formula <- "null"
  smiles <- "null"
  inchi<- "null"
  chemical_name<- "null"
  uri<- "null"
  theoretical_neutral_mass<- "null"
  adduct_ions <- paste(x$annotation,collapse = sep_val)
  reliability <- "null"
  best_id_confidence_measure <- "null"
  abundances <- x[,..idx_abundance]
  cnames <- colnames(abundances)
  quant <- apply(abundances,2,sum)
  unique_ms2 <- unique(x$opt_global_mgf_index_ms2)
  ms2_str <- NULL
  if(unique_ms2=="null"){
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

cnames <- colnames(sml_table)
sml_table <- cbind(rep("SML",nrow(sml_table)),1:nrow(sml_table),sml_table)
colnames(sml_table) <- c("SMH","SML_ID",cnames)



fwrite(sml_table,PATH_OUTPUT,sep="\t",append = TRUE,col.names = TRUE)
f <- file(PATH_OUTPUT,"a")
writeLines(c(""),f)
close(f)
fwrite(smf_table,PATH_OUTPUT,sep="\t",append = TRUE,col.names = TRUE)
