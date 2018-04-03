library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE89134/"

expr_id   <- "syn12063488"
anno_id   <- "syn12063555"
upload_id <- "syn12063036"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()



anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 31, nrows = 41) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[c(7,9), 1:5] %>% 
    transpose_df("attribute", "sample") %>% 
    set_colnames(c("sample", "cell_type", "patient")) %>% 
    mutate(patient = str_match(patient, "donor id: ([:alnum:]+)$")[,2]) %>% 
    mutate(cell_type = str_replace_all(as.character(cell_type), "ï", "i")) %>% 
    mutate(sample = str_replace_all(sample, "ï", "i")) 


expr_df <- expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("Hugo" = Gene) %>% 
    .[,c(1,7,13,2,8)] %>% 
    set_colnames(c("Hugo", anno_df$sample))

activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE89134/create_processed_tables.R")
)

write_tsv(expr_df, "counts.tsv")
upload_file_to_synapse("counts.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

