###WE load all the peaktables eventually
library(onlineLCMSaligner)

PATH_PEAKTABLES <- file.path("inst/peaktables")

all_peaktables <- list.files(PATH_PEAKTABLES,full.name = TRUE)

###Unicity of the alignement



###We just tes tthe lainger quickly

lam <- onlineLCMSaligner:::LCMSAlignerModelFromDirectory("inst/peaktables",path_model="inst/ex_model.rds",output="inst/blocks",save_interval=2,
                              num_file=2,num_peaks=100,rt=0.02,col_int="intensity",ppm=10,dmz=0.005)


lamt <- onlineLCMSaligner:::LCMSAlignerModelFromDirectory("inst/peaktables",path_model="inst/ex_model_nn.rds",output="inst/blocks_nn",save_interval=2,
                                                         num_file=2,num_peaks=100,rt=0.02,col_int="intensity",ppm=10,dmz=0.005,algorith="nn")



###We compare the two values eventually
table(lam@peaks$num)
table(lamt@peaks$num)


lam <- readRDS("inst/ex_model.rds")
pt <- read.table("inst/peaktables/peaktable_5_1553396fc0f446fd784d6eccacd24205.csv",sep=",")


any(is.na(pt$V2))
