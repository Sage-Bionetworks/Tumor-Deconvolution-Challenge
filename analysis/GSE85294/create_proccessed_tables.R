library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE85294/"

expr_id   <- "syn12063795"
anno_id   <- "syn12063794"
upload_id <- "syn12063039"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()



anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 36, nrows = 20) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[c(9,15),] %>% 
    transpose_df("attribute", "sample") %>% 
    select(3,2,1) %>% 
    set_colnames(c("sample", "cell_type", "replicate")) %>% 
    mutate(cell_type = str_replace_all(as.character(cell_type), "cell type: ", "")) %>% 
    mutate(cell_type = str_match(cell_type, "^([^\\(]*)")[,2]) %>% 
    mutate(cell_type = str_sub(cell_type, end = -2)) %>% 
    mutate(replicate = str_sub(as.character(replicate), start = -1))

expr_df <- expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    separate(Feature, into = c("Ensembl", "Hugo"), sep = "_")

activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE85294/create_processed_tables.R")
)

write_tsv(expr_df, "counts.tsv")
upload_file_to_synapse("counts.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

