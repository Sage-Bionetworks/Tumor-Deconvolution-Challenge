library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)

upload_id     <- "syn17026017"
gt_upload_id  <- "syn17026016"
sdy_id        <- "SDY144"
lt_id         <- "syn13363369"

source("../../scripts/utils.R")
synLogin()

cell.pheno.cols <- c("BASOPHIL_percent", "EOSINOPHIL_percent", "LYMPHOCYTE_percent",
                     "MONOCYTE_percent", "NEUTROPHIL_percent")

## Omit any NAs from lab tests ground truth and ensure that percents sum (nearly) to 100
ground_truth_lab_tests_df <- lt_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession %in% sdy_id) %>% 
    dplyr::rename(sample = subject_accession) %>%
    dplyr::select("sample", "age", "gender", "race", !!cell.pheno.cols) %>%
    na.omit() %>%
    .[rowSums(.[, c("BASOPHIL_percent", "EOSINOPHIL_percent", "LYMPHOCYTE_percent",
                    "MONOCYTE_percent", "NEUTROPHIL_percent")]) %in% c(100,101), ] %>%
    mutate(sample = paste0(sample, ".144"))		    
    

connection <- ImmuneSpaceR::CreateConnection(sdy_id)

expr_obj <- connection$getGEMatrix("SDY144_Other_TIV_Geo")
## expr_obj <- connection$getGEMatrix("SDY144_TIV2011_geo")

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

samples_in_common <- intersect(expr_df$sample, ground_truth_lab_tests_df$sample)

expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    transpose_df("sample", "Hugo")

ground_truth_df <- ground_truth_lab_tests_df %>% 
    filter(sample %in% samples_in_common) %>% 
    dplyr::select(sample, !!cell.pheno.cols)

anno_df <- ground_truth_lab_tests_df %>% 
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

upload_file_to_synapse("expression_illumina.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("annotation.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("ground_truth.tsv", gt_upload_id, activity_obj = activity_obj)

