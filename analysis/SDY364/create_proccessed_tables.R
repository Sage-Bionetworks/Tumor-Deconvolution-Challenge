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


connection <- ImmuneSpaceR::CreateConnection("SDY364")
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
    filter(study_accession == "SDY364") %>% 
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
    dplyr::filter(study_accession == "SDY364") %>% 
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



manifest_df1 <- tibble(
    path = c("expression_log.tsv",
             "expression_linear.tsv"),
    parent = upload_id,
    used = str_c(blood_count_id, expr_id, sep = ";"),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY364/create_processed_tables.R",
    activityName = "Process files from 10kimmunome.",
    file_type = c(rep("expression", 2)),
    expression_type = c(rep("microarray", 2)),
    microarray_type = c(rep("unknown", 2)),
    expression_space = c("log2", "linear")
)

manifest_df2 <- tibble(
    path = c("annotation.tsv"),
    parent = upload_id,
    used = str_c(blood_count_id, expr_id, sep = ";"),
    executed = "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY364/create_processed_tables.R",
    activityName = "Process files from 10kimmunome.",
    file_type = "annotations",
    annotations = "age;race;gender"
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


write_tsv(manifest_df1, "manifest1.tsv")
write_tsv(manifest_df2, "manifest2.tsv")
write_tsv(manifest_df3, "manifest3.tsv")

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(annotation_df, "annotation.tsv")
write_tsv(fcs_df, "ground_truth2.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

syncToSynapse("manifest1.tsv")
syncToSynapse("manifest2.tsv")
syncToSynapse("manifest3.tsv")

