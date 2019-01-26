library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)

upload_id     <- "syn17026017"
gt_upload_id  <- "syn17026016"
sdy_id        <- "SDY144"
lt_id         <- "syn13363369"
dataset       <- "GSE52005"
script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE52005/create_processed_tables.R"
activity_name <- "process_files_from_ImmuneSpace"


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
    as_tibble() %>%
    filter(study_time_collected == 0) %>% 
    dplyr::select(expr_id, participant_id) %>% 
    dplyr::rename(sample = participant_id) 

expr_df <- expr_obj@assayData$exprs %>% 
    t %>% 
    matrix_to_df("expr_id") %>% 
    inner_join(translation_df) %>% 
    dplyr::select(-expr_id)

samples_in_common <- intersect(expr_df$sample, ground_truth_lab_tests_df$sample)

ground_truth_df <- ground_truth_lab_tests_df %>% 
    filter(sample %in% samples_in_common) %>% 
    dplyr::select(sample, !!cell.pheno.cols)

anno_df <- ground_truth_lab_tests_df %>% 
    filter(sample %in% samples_in_common) %>% 
    dplyr::select(sample, age, gender, race) %>% 
    distinct

expr_df2 <- expr_df %>%
    filter(sample %in% samples_in_common) %>%
    gather(key = "gene", value = "expr", -sample)

linear_expr_df <- expr_df2 %>% 
    spread(key = "sample", value = "expr")

log_expr_df <- expr_df2 %>% 
    mutate(expr = log2(expr +1)) %>% 
    spread(key = "sample", value = "expr")

remove(expr_obj, expr_df, expr_df2)

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")
write_tsv(anno_df, "annotation.tsv")


expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    used = lt_id, 
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Illumina BeadArray Reader 500X",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    used = lt_id, 
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "annotations",
    annotations = "age;race;gender"
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    used = lt_id, 
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "percent",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)


write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")
