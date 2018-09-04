#open as sudo

library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(affy)
library(annotate)
library(org.Hs.eg.db)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE65135/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


# brian array file
brain_array_id <- "syn16781610"
brain_array_file <- download_from_synapse(brain_array_id)
command1 <- str_c("R CMD INSTALL ", brain_array_file)
system(command1)

# cel tar
cel_id <- "syn16781603"
cel_tar <- download_from_synapse(cel_id)
file.copy(cel_tar, basename(cel_tar))
command2 <- str_c("tar -xvf ", basename(cel_tar))
system(command2)


# from https://cibersort.stanford.edu/download.php ------
library(hgu133plus2hsentrezgcdf)
Data<-ReadAffy(cdfname = "hgu133plus2hsentrezgcdf")
eset<-mas5(Data)
ID<-featureNames(eset)
ID2<-sub("_at","",ID)
GS <- as.matrix(getSYMBOL(ID2, 'org.Hs.eg'))
ematrix<-exprs(eset)
rows <- GS
cols =  c("GeneSymbol",colnames(ematrix))
ematrix <- cbind(rows,ematrix)
ematrix <- ematrix[which(ematrix[,1] != "NA"),] #remove NAs
ematrix <- ematrix[order(ematrix[,1]),] #sort by gene name 
ematrix <- rbind(cols, ematrix)
# -----

expr_df <- ematrix %>% 
    .[-1,] %>% 
    as_data_frame() %>% 
    set_colnames(str_colnames(.))

write.table(ematrix,file="NormalizedExpressionArray.customCDF.txt",sep="\t", col.names=F, row.names=F,quote=FALSE)

