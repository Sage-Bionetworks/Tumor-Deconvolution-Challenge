library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

upload_id     <- "syn17096664"
gt_upload_id  <- "syn17096663"

source("../../scripts/utils.R")
synLogin()

cytof_id <- "syn13363372"
expr_id  <- "syn13363367"

cytof_df <- cytof_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession == "SDY113") %>% 
    select(-study_accession) %>% 
    dplyr::rename(sample = subject_accession)

expr_df <- expr_id  %>% 
    download_from_synapse %>%
    fread %>% 
    data.frame %>% 
    filter(study_accession == "SDY113") %>% 
    dplyr::select(-c(study_accession, data_accession, age, gender, race)) %>% 
    dplyr::rename(sample = subject_accession)

samples_in_common <- intersect(cytof_df$sample, expr_df$sample)

expr_df <- expr_df %>%
    filter(sample %in% samples_in_common) 

ground_truth_df <- cytof_df %>%
    dplyr::select(-c(age, gender, race)) %>%
    filter(sample %in% samples_in_common)

annotation_df <- cytof_df %>%
    dplyr::select(sample, age, gender, race) %>%
    filter(sample %in% samples_in_common)

activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY113/create_processed_tables.R")
)

write_tsv(expr_df, "expression.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(annotation_df, "annotation.tsv")

upload_file_to_synapse("expression.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

