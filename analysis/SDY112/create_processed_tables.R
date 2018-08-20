library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)


gt_upload_id <- "syn15588640"
upload_id    <- "syn15588641"

sdy_ids <- c("SDY112")

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/SDY112/"

fc_wb_id      <- "syn13363370"
cytof_pbmc_id <- "syn13363372"
expr_wb_id    <- "syn13363367"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


expr_wb_df <- expr_wb_id %>%
    download_from_synapse %>%
    fread %>%
    data.frame %>%
    inset("expression_source", value = "whole_blood") %>% 
    filter(study_accession %in% sdy_ids) %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_flow_wb_df <- fc_wb_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession %in% sdy_ids) %>% 
    .[, colSums(is.na(.)) < nrow(.)] %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_cytof_pbmc_df <- cytof_pbmc_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession %in% sdy_ids) %>% 
    .[, colSums(is.na(.)) < nrow(.)] %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_df <- 
    left_join(
        ground_truth_cytof_pbmc_df, 
        ground_truth_flow_wb_df, 
        by = c("study_accession", "age", "gender", "race", "sample")) %>% 
    dplyr::rename(T_cells = T_cells.x) %>% 
    dplyr::rename(T_cells_wb = T_cells.y) %>% 
    filter(sample %in% expr_wb_df$sample)

anno_df <- ground_truth_df %>% 
    select(sample, age, race, gender)

ground_truth_df <- ground_truth_df %>% 
    select(-c(age, race, gender, study_accession))

remove(ground_truth_flow_wb_df, ground_truth_cytof_pbmc_df)

expr_df <- expr_wb_df %>% 
    filter(sample %in% ground_truth_df$sample) %>% 
    select(-c(study_accession, data_accession, age, gender, race)) %>% 
    as.data.frame %>% 
    column_to_rownames("sample") %>%
    data.matrix(rownames.force = NA) %>% 
    t %>% 
    as.data.frame %>%
    rownames_to_column("Hugo") %>%
    as_data_frame %>%
    group_by(Hugo) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    .[complete.cases(.), ]

remove(expr_wb_df)

activity_obj <- Activity(
    name = "create",
    description = "process 10Kimmunome files into usable tables",
    used = list(cytof_pbmc_id, fc_wb_id, expr_wb_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY112/create_processed_tables.R"))

write_tsv(expr_df, "expression_microarray.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

upload_file_to_synapse("expression_microarray.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)
