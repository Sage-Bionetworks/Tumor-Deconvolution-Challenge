library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(Biobase)
library(ImmuneSpaceR)

upload_id     <- "syn18349400"

expr_whole_blood_10K_id    <- "syn13363367"

source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY305/create_processed_tables.R"
dataset       <- "SDY305"
activity_name <- "process Immunespace and 10k immunomes files into usable tables"
used          <- str_c(expr_whole_blood_10K_id, sep = ";")

expr_wb_10K_df <- expr_whole_blood_10K_id %>%
    download_from_synapse %>%
    fread %>%
    data.frame %>% 
    filter(study_accession == dataset) %>% 
    select(-c(study_accession, data_accession, age, gender, race)) %>% 
    as_tibble() %>% 
    dplyr::rename(sample = subject_accession) %>% 
    distinct %>% 
    gather(key = "Hugo", value = "expr", -c("sample")) %>% 
    drop_na()


con   <- ImmuneSpaceR::CreateConnection(dataset)

expr_sets <- c(
    "SDY305_Other_IDTIV_Geo",
    "SDY305_Other_TIV_Geo"
)


expr_immunespace_df <- expr_sets %>%
    con$getGEMatrix() %>%
    con$mapSampleNames(EM = .) %>%
    Biobase::exprs() %>%
    matrix_to_df("Hugo") %>%
    gather(key = "sample", value = "expr",  -Hugo) %>%
    separate(sample, sep = "_", into = c("sample", "time")) %>%
    filter(time == "d0") %>%
    select(-time) %>%
    mutate(sample = str_remove_all(sample, "\\.[0-9]*"))

ground_truth_df <-
    con$getDataset("fcs_analyzed_result") %>%
    dplyr::as_tibble() %>%
    mutate(sample = str_sub(participant_id, end = -5)) %>%
    filter(study_time_collected == 0)


# samples_in_common <- intersect(expr_df$sample, ground_truth_df$sample)
# 
# 
# linear_expr_df <- expr_df %>%
#     filter(sample %in% samples_in_common) %>% 
#     spread(key = "sample", value = "expr") 
# 
# log_expr_df <- expr_df %>%
#     filter(sample %in% samples_in_common) %>% 
#     mutate(expr = log2(expr)) %>%
#     spread(key = "sample", value = "expr") 
# 
# ground_truth_df <- ground_truth_immune_space %>%
#     filter(sample %in% samples_in_common) %>% 
#     spread(key = "cell_type", value = "percent")
# 
# 
# 
# expression_manifest_df <- tibble(
#     path = c("expression_log.tsv", "expression_linear.tsv"),
#     parent = upload_id,
#     executed = script_url,
#     activityName = activity_name,
#     dataset = dataset,
#     used = used,
#     file_type = "expression",
#     expression_type = "microarray",
#     microarray_type = "Illumina Human HT-12 V3 BeadChip",
#     expression_space = c("log2", "linear")
# )
# 
# # is cell number the unit?
# ground_truth_manifest_df <- tibble(
#     path = "ground_truth.tsv",
#     parent = upload_id,
#     executed = script_url,
#     activityName = activity_name,
#     dataset = dataset,
#     used = used,
#     file_type = "ground truth",
#     unit = "cell number",
#     cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
# )
# 
# write_tsv(log_expr_df, "expression_log.tsv")
# write_tsv(linear_expr_df, "expression_linear.tsv")
# write_tsv(ground_truth_df, "ground_truth.tsv")
# 
# write_tsv(expression_manifest_df, "expression_manifest.tsv")
# write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")
# 
# syncToSynapse("expression_manifest.tsv")
# syncToSynapse("ground_truth_manifest.tsv")
