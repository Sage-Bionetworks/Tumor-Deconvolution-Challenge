library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)


gt_upload_id <- "syn13363373"
upload_id <- "syn13363374"

gt_pbmc_id   <- "syn13363371"
expr_pbmc_id <- "syn13363368"
anno_id      <- "syn13363390"

sdy_id <- "SDY404"
home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE59654/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

gt_pbmc_df <- gt_pbmc_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession == sdy_id) %>% 
    select(-study_accession) %>% 
    dplyr::rename(sample = subject_accession)

expr_pbmc_df <- expr_pbmc_id %>% 
    download_from_synapse %>%
    fread %>% 
    data.frame %>% 
    filter(study_accession == sdy_id) %>% 
    select(-study_accession) %>% 
    dplyr::rename(sample = subject_accession)

anno_pbmc_df <- expr_pbmc_df %>% 
    as_data_frame %>% 
    select(sample, data_accession, age, gender, race)


geo_anno_df <- anno_id %>% 
    create_df_from_synapse_id(unzip = T, skip = 34, nrow = 30) %>%
    dplyr::rename("title" = `!Sample_title`) %>% 
    .[c(1, 9, 10,11), ] %>% 
    mutate(title = c("data_accession", "day", "age_group", "phenotype")) %>% 
    transpose_df("title", "sample") %>% 
    filter(day == "blood draw date: 0") %>% 
    select(-c(day, sample)) %>% 
    mutate(age_group = str_remove_all(age_group, "age group: ")) %>% 
    mutate(phenotype = str_remove_all(phenotype, "phenotype: "))
    
anno_df <- geo_anno_df %>% 
    inner_join(anno_pbmc_df) %>% 
    inner_join(gt_pbmc_df[,1:4]) %>% 
    distinct %>% 
    select(-data_accession)

ground_truth_df <- gt_pbmc_df %>% 
    select(-c(age, race, gender)) %>% 
    group_by(sample) %>% 
    summarise_all(mean, na.rm = T) %>% 
    ungroup

expr_df <- expr_pbmc_df %>% 
    select(-c(data_accession, age, gender, race)) %>% 
    column_to_rownames("sample") %>% 
    as.matrix %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column("Hugo") %>% 
    as_data_frame %>% 
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup

activity_obj <- Activity(
    name = "create",
    description = "process GEO files into usable tables",
    used = list(gt_pbmc_id, expr_pbmc_id, anno_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/create_processed_tables.R")
)

write_tsv(expr_df, "expression_microarray.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(anno_df, "annotation.tsv")

upload_file_to_synapse("expression_microarray.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
