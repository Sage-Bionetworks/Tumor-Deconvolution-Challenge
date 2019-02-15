library(synapser)
library(data.table)
library(magrittr)
library(tidyverse)
library(Biobase)

sdy_ids <- c("SDY311", "SDY312", "SDY314", "SDY315")

fc_wb_id      <- "syn13363370"
expr_wb_id    <- "syn13363367"


source("../../scripts/utils.R")
synLogin()


# connection <- ImmuneSpaceR::CreateConnection("SDY311")
connection <- ImmuneSpaceR::CreateConnection("SDY312")
# connection <- ImmuneSpaceR::CreateConnection("SDY314")
# connection <- ImmuneSpaceR::CreateConnection("SDY315")

connection$listDatasets()
x <- connection$getGEMatrix(c("SDY312_Other_GroupA", "SDY312_Other_GroupB", "SDY312_Other_GroupC"))
x2 <- exprs(connection$mapSampleNames(EM = x))


str_sub(colnames(x2), end = -8) %in% ground_truth_flow_wb_df$sample
ground_truth_flow_wb_df$sample %in% str_sub(colnames(x2), end = -8)





expr_wb_df <- expr_wb_id %>%
    download_from_synapse %>%
    fread %>%
    data.frame %>%
    inset("expression_source", value = "whole_blood") %>%
    filter(study_accession %in% sdy_ids)


ground_truth_flow_wb_df <- fc_wb_id %>%
    create_df_from_synapse_id %>%
    filter(study_accession %in% sdy_ids) %>%
    .[, colSums(is.na(.)) < nrow(.)] %>%
    dplyr::rename(sample = subject_accession)


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
