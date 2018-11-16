library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)

upload_id     <- "syn17026017"
gt_upload_id  <- "syn17026016"
sdy_id        <- "SDY144"

source("../../scripts/utils.R")
synLogin()



connection <- ImmuneSpaceR::CreateConnection(sdy_id)
flow_df    <- connection$getDataset("fcs_analyzed_result") %>% 
    as_data_frame() %>% 
    filter(study_time_collected == 0) %>% 
    filter(cell_number_unit == "cells/uL") %>% 
    dplyr::select("participant_id", "age_reported", "gender", "race", "population_cell_number", 
           "cell_number_unit", "population_definition_reported", "population_name_reported") %>% 
    set_colnames(c("sample", "age", "gender", "race", "population_cell_number", 
                   "cell_number_unit", "population_definition_reported", "population_name_reported"))

expr_obj <- connection$getGEMatrix("SDY144_TIV2011_geo")

translation_df <- expr_obj@phenoData@data %>% 
    rownames_to_column("expr_id") %>% 
    as_data_frame() %>%
    filter(study_time_collected == 0) %>% 
    dplyr::select(expr_id, participant_id) %>% 
    dplyr::rename(sample = participant_id)

expr_df <- expr_obj@assayData$exprs %>% 
    t %>% 
    matrix_to_df("expr_id") %>% 
    inner_join(translation_df) %>% 
    dplyr::select(-expr_id)

samples_in_common <- intersect(expr_df$sample, flow_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    transpose_df("sample", "Hugo")

ground_truth_df <- flow_df %>% 
    filter(sample %in% samples_in_common) %>% 
    dplyr::select(sample, population_cell_number, population_name_reported) %>% 
    spread(key = population_name_reported, value = population_cell_number)

anno_df <- flow_df %>% 
    filter(sample %in% samples_in_common) %>% 
    dplyr::select(sample, age, gender, race) %>% 
    distinct


activity_obj <- Activity(
    name = "create",
    description = "process GEO data into usable tables",
    used = list(),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE52005/create_processed_tables.R")
)

write_tsv(expr_df, "expression_illumina.tsv")
write_tsv(anno_df, "annotation.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(flow_df, "flow.tsv")

upload_file_to_synapse("expression_illumina.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)
upload_file_to_synapse("flow.tsv", upload_id, activity_obj = activity_obj)

