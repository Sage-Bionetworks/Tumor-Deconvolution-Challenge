library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE89225/"

ill_expr_id <- "syn12063041"
ion_expr_id <- "syn12063040"
ill_anno_id <- "syn12063043"
ion_anno_id <- "syn12063042"
hugo_id     <- "syn11536071"
upload_id   <- "syn12063034"


hugo_df <-  create_df_from_synapse_id(hugo_id)

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()



ion_anno_df <- ion_anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 36, nrows = 41) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[c(7,9,36),] %>% 
    transpose_df("attribute", "sample") %>% 
    set_colnames(c("sample", "cell_type", "source", "platform")) %>% 
    mutate(patient = str_match(sample, "_([:alnum:]+)$")[,2])

ill_anno_df <- ill_anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 36, nrows = 41) %>% 
    dplyr::rename("attribute" = `!Sample_title`) %>% 
    .[c(7,9,37),] %>% 
    transpose_df("attribute", "sample") %>% 
    set_colnames(c("sample", "cell_type", "source", "platform")) %>% 
    mutate(patient = str_match(sample, "_([:alnum:]+)$")[,2])

anno_df <- bind_rows(ion_anno_df, ill_anno_df)



ion_expr_df <- ion_expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("Hugo" = V1)

ill_expr_df <- ill_expr_id %>% 
    create_df_from_synapse_id(unzip = T) %>% 
    dplyr::rename("ensembl_gene_id" = V1) %>% 
    left_join(hugo_df) %>%
    dplyr::rename("Hugo" = hgnc_symbol) %>% 
    select(-ensembl_gene_id) %>% 
    filter(!Hugo == "") %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup

expr_df <- inner_join(ion_expr_df, ill_expr_df)
    
            







files <- list.files() %>% 
    keep(str_detect(., "RPKMforgenes.txt$"))

sample_names <- files %>% 
    str_split("_") %>% 
    map_chr(extract2, 1)

count_dfs <- files %>% 
    map(fread) %>% 
    map(as_data_frame) %>% 
    map(select, c(1,4)) %>% 
    map(set_colnames, c("Hugo", "Count")) 

genes <- count_dfs[[1]]$Hugo

count_df <- count_dfs %>% 
    map(select, 2) %>% 
    bind_cols %>% 
    set_colnames(sample_names) %>% 
    inset("Hugo", value = genes) %>% 
    select(Hugo, everything())



activity_obj <- Activity(
    name = "create",
    description = "process GEO file into usable table",
    used = list(expr_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE89225/create_processed_tables.R")
)

write_tsv(count_df, "counts.tsv")
upload_file_to_synapse("counts.tsv", upload_id, activity_obj = activity_obj)

write_tsv(anno_df, "annotation.tsv")
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)

