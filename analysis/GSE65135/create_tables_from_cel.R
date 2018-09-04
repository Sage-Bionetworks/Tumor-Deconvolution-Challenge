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

upload_id <- "syn16781521"
brain_array_id <- "syn16781610"
cel_id <- "syn16781603"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


# brian array file
brain_array_file <- download_from_synapse(brain_array_id)
command1 <- str_c("R CMD INSTALL ", brain_array_file)
system(command1)

# cel tar
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
    set_rownames(.[,1]) %>% 
    .[,-1] %>% 
    set_colnames(str_match(colnames(.), "(GSM[0-9]+)_[:print:]+")[,2]) %>% 
    data.frame() %>% 
    rownames_to_column("Hugo") %>% 
    as_data_frame() %>% 
    gather(key = "sample", value = "expr", -Hugo) %>% 
    mutate(expr = as.numeric(expr))

expr_df %>% 
    spread(key = "sample", value = "expr") %>% 
    write_tsv("expression_affy.tsv")


expr_df %>% 
    mutate(expr = log2(expr)) %>% 
    spread(key = "sample", value = "expr") %>% 
    write_tsv("log_expression_affy.tsv")

activity_obj <- Activity(
    name = "create",
    description = "process GEO CEL data into usable tables",
    used = list(brain_array_id, cel_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65135/create_tables_from_cel.R")
)

upload_file_to_synapse("log_expression_affy.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("expression_affy.tsv", upload_id, activity_obj = activity_obj)