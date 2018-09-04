library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)


gt_upload_id <- "syn14566979"
upload_id    <- "syn14567262"

fc_pbmc_id    <- "syn13363371"
fc_wb_id      <- "syn13363370"
expr_wb_id    <- "syn13363367"

sdy_ids  <- c("SDY212")
home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE41080/"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()

ground_truth_pbmc_df <- fc_pbmc_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession %in% sdy_ids) %>% 
    .[, colSums(is.na(.)) < nrow(.)] %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_wb_df <- fc_wb_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession %in% sdy_ids) %>% 
    .[, colSums(is.na(.)) < nrow(.)] %>% 
    dplyr::rename(sample = subject_accession)

expr_wb_df <- expr_wb_id %>%
    download_from_synapse %>%
    fread %>%
    data.frame %>%
    inset("expression_source", value = "whole_blood") %>% 
    filter(study_accession %in% sdy_ids) 
    
expr_df <- expr_wb_df %>% 
    dplyr::rename(sample = subject_accession) %>% 
    as_data_frame
    
remove(expr_pbmc_df, expr_wb_df)
 
anno_expr_df <- select(
    expr_df, data_accession, sample, age, gender, race, study_accession, expression_source)

anno_gt_df <- select(ground_truth_pbmc_df, sample, age, gender, race, study_accession)

anno_df <- 
    inner_join(anno_expr_df, anno_gt_df) %>% 
    select(sample, age, gender, race) %>% 
    distinct

expr_df <- expr_df %>% 
    select(-c(study_accession, data_accession, age, gender, race, expression_source)) %>% 
    as.data.frame %>% 
    column_to_rownames("sample") %>%
    data.matrix(rownames.force = NA) %>% 
    t %>% 
    as.data.frame %>%
    rownames_to_column("Hugo") %>%
    as_data_frame %>%
    group_by(Hugo) %>% 
    summarise_all(max) %>% 
    ungroup %>% 
    select(c("Hugo", anno_df$sample))

ground_truth_pbmc_df <- ground_truth_pbmc_df %>% 
    select(- c(study_accession, age, race, gender)) %>%
    filter(sample %in% anno_df$sample) 

ground_truth_wb_df <- ground_truth_wb_df %>% 
    select(- c(study_accession, age, race, gender)) %>%
    filter(sample %in% anno_df$sample) 

activity_obj <- Activity(
    name = "create",
    description = "process 10Kimmunome files into usable tables",
    used = list(fc_pbmc_id, fc_wb_id, expr_wb_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE41080/create_processed_tables.R")
)

write_tsv(expr_df, "expression_microarray.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_pbmc_df, "ground_truth_pbmc.tsv")
write_tsv(ground_truth_wb_df, "ground_truth_wb.tsv")


upload_file_to_synapse("expression_microarray.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth_pbmc.tsv", gt_upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth_wb.tsv", gt_upload_id, activity_obj = activity_obj)

