library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)
library(synapserutils)

gt_upload_id  <- "syn17971750"

## preprocessed_folder_upload_id, gt_upload_id, and dataset defined in setup.R
source("setup.R")
source("../../scripts/utils.R")
synLogin()

## 10KImmunomes.Lab Tests_ Blood Count
lab_tests_10k_id <- "syn13363369"

## 10KImmunomes.Gene Expression_ Whole Blood
expr_wb_10k_id  <- "syn13363367"

synLogin()

dataset <- "SDY364"
github.path   <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path  <- paste0(github.path, dataset)
script_url    <- paste0(dataset.path, "/create_processed_tables.R")
used          <- NA
activity_name <- "process GEO data into usable tables"

connection <- ImmuneSpaceR::CreateConnection(dataset)
fcs_df <- 
    connection$getDataset("fcs_analyzed_result") %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(study_time_collected == 0) %>% 
    mutate(sample = str_sub(participant_id, end = -5))

if(any(fcs_df$cell_number_unit != "ul")) {
  stop("Was expecting all ground truth units to be ul.\n")
}

expr_and_anno_df <- expr_wb_10k_id  %>% 
    download_from_synapse %>%
    fread %>% 
    data.frame %>% 
    dplyr::filter(study_accession == dataset) %>%
    dplyr::rename(sample = subject_accession) 

expr_df <- expr_and_anno_df %>%
    dplyr::select(-c(study_accession, data_accession, age, gender, race)) %>% 
    as_tibble() %>% 
    gather(key = "Hugo", value = "expr", - sample) %>% 
    drop_na()

## Output annotation, ground truth, and expression for flow data
postfix <- "_flow"

annotation_df <- expr_and_anno_df %>%
    dplyr::select(c(sample, age, gender, race)) %>%
    unique()

ground_truth_df <- fcs_df %>% 
    select(sample, population_name_reported, population_cell_number) %>%
    spread(key = "population_name_reported", value = "population_cell_number")    

samples_in_common <- reduce(list(expr_df$sample, ground_truth_df$sample, annotation_df$sample), intersect)

linear_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

log_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    mutate(expr = log2(expr +1)) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common)

annotation_df <- annotation_df %>%
    filter(sample %in% samples_in_common)

expression_manifest_df <- tibble(
    path = c(paste0("expression_log", postfix, ".tsv"), paste0("expression_linear", postfix, ".tsv")),
    parent = preprocessed_folder_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Illumina BeadArray",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = paste0("annotation", postfix, ".tsv"),
    parent = preprocessed_folder_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = paste0("ground_truth", postfix, ".tsv"),
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "ul",
    cell_types = c(str_c(colnames(ground_truth_df)[-1], collapse = ";"))
)

write_tsv(log_expr_df, paste0("expression_log", postfix, ".tsv"))
write_tsv(linear_expr_df, paste0("expression_linear", postfix, ".tsv"))
write_tsv(annotation_df, paste0("annotation", postfix, ".tsv"))
write_tsv(ground_truth_df, paste0("ground_truth", postfix, ".tsv"))

write_tsv(expression_manifest_df, paste0("expression_manifest", postfix, ".tsv"))
write_tsv(annotation_manifest_df, paste0("annotation_manifest", postfix, ".tsv"))
write_tsv(ground_truth_manifest_df, paste0("ground_truth_manifest", postfix, ".tsv"))

syncToSynapse(paste0("expression_manifest", postfix, ".tsv"))
syncToSynapse(paste0("annotation_manifest", postfix, ".tsv"))
syncToSynapse(paste0("ground_truth_manifest", postfix, ".tsv"))


## Output annotation, ground truth, and expression for lab tests data
postfix <- "_labs"

annotation_df <- expr_and_anno_df %>%
    dplyr::select(c(sample, age, gender, race)) %>%
    unique()

blood_count_df <- lab_tests_10k_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession == dataset) %>% 
    select(-c(study_accession, WBC_K_per_uL)) %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_df <- blood_count_df %>% 
    dplyr::select(sample, BASOPHIL_percent, EOSINOPHIL_percent, LYMPHOCYTE_percent, MONOCYTE_percent, NEUTROPHIL_percent) %>%
    na.omit()

unity <- 100
eps <- unity * 0.01

flag <- (abs(ground_truth_df %>% df_to_matrix("sample") %>% rowSums() - unity) > eps)
if(any(flag)) {
  cat("The following rows in the ground truth do not sum to one:\n")
  print(ground_truth_df[flag,,drop=F])
}

## Convert ground truth to fraction
ground_truth_df <- ground_truth_df %>%
  df_to_matrix("sample") %>%
  divide_by(rowSums(.)) %>%
  matrix_to_df("sample")

samples_in_common <- reduce(list(expr_df$sample, ground_truth_df$sample, annotation_df$sample), intersect)

linear_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

log_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    mutate(expr = log2(expr +1)) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common)

annotation_df <- annotation_df %>%
    filter(sample %in% samples_in_common)

expression_manifest_df <- tibble(
    path = c(paste0("expression_log", postfix, ".tsv"), paste0("expression_linear", postfix, ".tsv")),
    parent = preprocessed_folder_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "Illumina BeadArray",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = paste0("annotation", postfix, ".tsv"),
    parent = preprocessed_folder_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

ground_truth_manifest_df <- tibble(
    path = paste0("ground_truth", postfix, ".tsv"),
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = c(str_c(colnames(ground_truth_df)[-1], collapse = ";"))
)

write_tsv(log_expr_df, paste0("expression_log", postfix, ".tsv"))
write_tsv(linear_expr_df, paste0("expression_linear", postfix, ".tsv"))
write_tsv(annotation_df, paste0("annotation", postfix, ".tsv"))
write_tsv(ground_truth_df, paste0("ground_truth", postfix, ".tsv"))

write_tsv(expression_manifest_df, paste0("expression_manifest", postfix, ".tsv"))
write_tsv(annotation_manifest_df, paste0("annotation_manifest", postfix, ".tsv"))
write_tsv(ground_truth_manifest_df, paste0("ground_truth_manifest", postfix, ".tsv"))

syncToSynapse(paste0("expression_manifest", postfix, ".tsv"))
syncToSynapse(paste0("annotation_manifest", postfix, ".tsv"))
syncToSynapse(paste0("ground_truth_manifest", postfix, ".tsv"))

