library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE103322/"

expr_id <- "syn11990085"
anno_id <- "syn11990149"

upload_id  <- "syn11990084"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

expr_raw_df <- create_df_from_synapse_id(expr_id, unzip = T)

expr_df <- expr_raw_df %>% 
    .[-c(1:5), ] %>% 
    rename("Hugo" = V1) %>% 
    mutate(Hugo = str_replace_all(Hugo, "'", "")) %>% 
    arrange(Hugo) %>% 
    .[,order(colnames(.))] %>% 
    select(Hugo, everything())

cell_type_df <- expr_raw_df %>% 
    .[5, -1] %>% 
    data.frame %>% 
    set_rownames("cell_type") %>% 
    as.matrix %>% 
    t %>% 
    matrix_to_df("sample") %>% 
    arrange(sample) %>% 
    mutate(cell_type = as.character(cell_type)) %>% 
    mutate(cell_type = ifelse(
        cell_type == "0", 
        "tumor", 
        ifelse(cell_type == "-Fibroblast",
               "Fibroblast",
               cell_type)))

anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 27, nrows = 35) %>% 
    .[10, -1] %>% 
    data.frame %>% 
    set_rownames("cell_source") %>% 
    as.matrix %>% 
    t %>% 
    matrix_to_df("sample") %>% 
    mutate(cell_source = str_sub(cell_source, start = 13)) %>% 
    left_join(cell_type_df) %>% 
    .[complete.cases(.),] %>% 
    filter(sample %in% colnames(expr_df))

expr_df <- expr_df[,colnames(expr_df) %in% c("Hugo", anno_df$sample)]

activity_obj <- Activity(
    name = "create",
    description = "process GEO files into usable tables",
    used = list(anno_id, expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE103322/create_processed_tables.R")
)

write_tsv(expr_df, "expression.tsv")
write_tsv(anno_df, "annotation.tsv")

upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)