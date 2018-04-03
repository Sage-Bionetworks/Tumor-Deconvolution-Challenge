library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE80306/"

expr_id   <- "syn12063785"
anno_id   <- "syn12063786"
upload_id <- "syn12063038"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()



anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 26, nrows = 33) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[9:10, ] %>% 
    mutate(attribute = c("patient", "cell_type")) %>% 
    transpose_df("attribute", "sample") %>% 
    set_colnames(c("sample", "patient", "cell_type")) %>% 
    mutate(cell_type = str_replace_all(as.character(cell_type), "cell type: ", "")) %>% 
    mutate(patient = str_replace_all(as.character(patient), "donor id: ", ""))


expr_df <- expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("Hugo" = gene) %>% 
    set_colnames(str_replace_all(colnames(.), "Hs", "HS"))

activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE80306/create_processed_tables.R")
)

write_tsv(expr_df, "counts.tsv")
upload_file_to_synapse("counts.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

