library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)
library(synapserutils)

upload_id     <- "syn17971758"
gt_upload_id  <- "syn17971750"

source("../../scripts/utils.R")
synLogin()

blood_count_id <- "syn13363369"
expr_id  <- "syn13363367"

dataset <- "SDY364"
script_url <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY364/create_processed_tables.R"
activity_name <- "Process files from 10kimmunome."
used <- str_c(blood_count_id, expr_id, sep = ";")



connection <- ImmuneSpaceR::CreateConnection(dataset)
fcs_df <- 
    connection$getDataset("fcs_analyzed_result") %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(study_time_collected == 0) %>% 
    mutate(sample = str_sub(participant_id, end = -5)) %>% 
    select(-c(
        participant_id,
        age_reported, 
        gender,
        race, 
        cohort, 
        study_time_collected, 
        study_time_collected_unit,
        study_time_collected_unit, 
        base_parent_population)) %>% 
    select(sample, everything())
        
    


blood_count_df <- blood_count_id %>% 
    create_df_from_synapse_id %>% 
    filter(study_accession == dataset) %>% 
    select(-c(study_accession, WBC_K_per_uL)) %>% 
    dplyr::rename(sample = subject_accession)



columns_with_data <- blood_count_df %>% 
    is.na() %>% 
    colSums() %>% 
    .[. < 23] %>% 
    names

blood_count_df <- select(blood_count_df, columns_with_data)

ground_truth_df <- blood_count_df %>% 
    dplyr::select(-c(age, race, gender)) %>% 
    gather(key = "celltype", value = "percent", -sample) %>% 
    mutate(fraction = percent / 100) %>% 
    mutate(celltype = str_remove_all(celltype, "_percent")) %>% 
    select(-percent) %>% 
    drop_na() %>% 
    spread(key = "celltype", value = "fraction")
    

annotation_df <- blood_count_df %>% 
    select(sample, age, race, gender) 




expr_df <- expr_id  %>% 
    download_from_synapse %>%
    fread %>% 
    data.frame %>% 
    dplyr::filter(study_accession == dataset) %>% 
    dplyr::select(-c(study_accession, data_accession, age, gender, race)) %>% 
    dplyr::rename(sample = subject_accession) %>% 
    as_tibble() %>% 
    gather(key = "Hugo", value = "expr", - sample) %>% 
    drop_na()

samples_in_common <- intersect(ground_truth_df$sample, expr_df$sample)


linear_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    spread(key = "sample", value = "expr")

log_expr_df <- expr_df %>% 
    filter(sample %in% samples_in_common) %>% 
    mutate(expr = log2(expr +1)) %>% 
    spread(key = "sample", value = "expr")
    


remove(expr_df)


ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common)

annotation_df <- annotation_df %>%
    filter(sample %in% samples_in_common)

fcs_df <- fcs_df %>% 
    filter(sample %in% samples_in_common)



expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "microarray", 
    microarray_type = "unknown",
    expression_space = c("log2", "linear")
)

annotation_manifest_df <- tibble(
    path = "annotation.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "annotations",
    annotations = str_c(colnames(annotation_df)[-1], collapse = ";")
)

manifest_df3 <- tibble(
    path = c("ground_truth.tsv", "ground_truth2.tsv"),
    parent = gt_upload_id,
    used = str_c(blood_count_id, expr_id, sep = ";"),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY364/create_processed_tables.R",
    activityName = "Process files from 10kimmunome.",
    file_type = "ground truth",
    unit = c("fraction", "ul"),
    cell_types = c("BASOPHIL;EOSINOPHIL;LYMPHOCYTE;MONOCYTE;NEUTROPHIL", "various")
)

ground_truth_manifest_df <- tibble(
    path = c("ground_truth.tsv", "ground_truth2.tsv"),
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = c("fraction", "ul"),
    cell_types = c(str_c(colnames(ground_truth_df)[-1], collapse = ";"), "various")
)

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(annotation_df, "annotation.tsv")
write_tsv(fcs_df, "ground_truth2.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(annotation_manifest_df, "annotation_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("annotation_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")

