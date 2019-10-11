library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)


sdy_ids <- c(
    "SDY296",
    "SDY301", 
    "SDY387", 
    "SDY522", 
    "SDY984"
    )

sdy_ids <- "SDY80"

fc_wb_id      <- "syn13363370"
expr_wb_id    <- "syn13363367"


source("../../scripts/utils.R")
synLogin()




expr_wb_df <- expr_wb_id %>%
    download_from_synapse %>%
    fread %>%
    data.frame %>%
    inset("expression_source", value = "whole_blood") %>%
    filter(study_accession %in% sdy_ids)


ground_truth_flow_wb_df <- fc_wb_id %>%
    create_df_from_synapse_id %>%
    filter(study_accession %in% sdy_ids) %>%
    .[, colSums(is.na(.)) < nrow(.)]


ground_truth_cytof_pbmc_df <- "syn13363372" %>%
    create_df_from_synapse_id %>%
    filter(study_accession %in% sdy_ids)
# 
# 
# 
# 
# 
# 
# 
# cell_type_df <- bind_rows(cytof_pbmc_df, flow_pbmc_df, flow_whole_blood_df)
# 
# cell_type_anno_df <- select(cell_type_df , study_accession, age, gender, race, subject_accession)
# 
# cell_type_df <- select(cell_type_df, -c(study_accession, age, gender, race))
# 
# 
# expr_pbmc_df <- "10KImmunomes.Gene Expression_ PBMC.2018-06-25.csv" %>% 
#     fread %>%
#     data.frame %>% 
#     filter(study_accession %in% studies)
# 
# anno_pbmc_df <- expr_pbmc_df %>% 
#     as_data_frame %>% 
#     select(subject_accession, data_accession, age, gender, race, study_accession)
# 
# expr_pbmc_df <- expr_pbmc_df %>% 
#     select(-c(data_accession, age, gender, race, study_accession)) %>% 
#     column_to_rownames("subject_accession") %>% 
#     as.matrix %>% 
#     t %>% 
#     as.data.frame %>% 
#     rownames_to_column("Hugo") %>% 
#     as_data_frame
# 
# expr_wb_df <- "10KImmunomes.Gene Expression_ Whole Blood.2018-06-25.csv" %>% 
#     fread %>%
#     data.frame %>% 
#     filter(study_accession %in% studies)
# 
# anno_wb_df <- expr_wb_df %>% 
#     as_data_frame %>% 
#     select(subject_accession, data_accession, age, gender, race, study_accession)
# 
# expr_wb_df <- expr_wb_df %>% 
#     select(-c(data_accession, age, gender, race, study_accession)) %>% 
#     column_to_rownames("subject_accession") %>% 
#     as.matrix %>% 
#     t %>% 
#     as.data.frame %>% 
#     rownames_to_column("Hugo") %>% 
#     as_data_frame
# 
# expr_anno_df <- bind_rows(anno_pbmc_df, anno_wb_df)
# 
# 
# 
