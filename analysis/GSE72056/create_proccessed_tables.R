library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE72056/"

expr_id   <- "syn11808157"
anno_id   <- "syn12117735"
upload_id <- "syn11808146"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

expr_df <- create_df_from_synapse_id(expr_id) %>% 
    set_colnames(str_replace_all(colnames(.), "-", "_")) %>% 
    set_colnames(str_replace_all(colnames(.), "\\.", "_")) 

# malignant(1=no,2=yes,0=unresolved)
# "non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"
expr_sample_df <- expr_df %>% 
    .[1:3,] %>% 
    transpose_df("Cell", "sample") %>% 
    set_colnames(c("sample", "patient", "cancer_status", "cell_type")) %>% 
    mutate(patient = str_c("CY", patient)) %>% 
    mutate(cancer_status = as.character(cancer_status)) %>% 
    mutate(cancer_status = ifelse(cancer_status == "1", "no", cancer_status)) %>%
    mutate(cancer_status = ifelse(cancer_status == "2", "yes", cancer_status)) %>%
    mutate(cancer_status = ifelse(cancer_status == "0", "unresolved", cancer_status)) %>%
    mutate(cell_type = as.character(cell_type)) %>% 
    mutate(cell_type = ifelse(cell_type == "0", "melanoma/unresolved", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "1", "T", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "2", "B", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "3", "Macro.", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "4", "Endo.", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "5", "CAF", cell_type)) %>%
    mutate(cell_type = ifelse(cell_type == "6", "NK", cell_type))

expr_sample_df1 <- expr_sample_df %>% 
    filter(cancer_status == "no") %>% 
    filter(cell_type != "melanoma/unresolved") %>% 
    select(-cancer_status)

expr_sample_df2 <- expr_sample_df %>% 
    filter(cancer_status == "yes") %>% 
    filter(cell_type == "melanoma/unresolved") %>% 
    mutate(cell_type = "melanoma") %>% 
    select(-cancer_status)

expr_sample_df <- bind_rows(expr_sample_df1, expr_sample_df2)

anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 20, nrows = 30) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[3,] %>% 
    transpose_df("attribute", "sample") %>% 
    set_colnames(c("sample", "batch")) %>% 
    mutate(sample = str_replace_all(sample, "-", "_")) %>% 
    mutate(sample = str_replace_all(sample, "\\.", "_")) %>% 
    inner_join(expr_sample_df) %>% 
    arrange(sample)

expr_df <- expr_df %>% 
    .[-c(1:3),] %>% 
    rename("Hugo" = Cell) %>% 
    .[,c("Hugo", anno_df$sample)]

activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE72056/create_processed_tables.R")
)

write_tsv(expr_df, "expression.tsv")
upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

