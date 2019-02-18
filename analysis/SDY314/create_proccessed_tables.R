library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(Biobase)
library(ImmuneSpaceR)

upload_id     <- "syn18347744"

import_expr_id    <- "syn18347746"
import_trans_id   <- "syn18347747"

source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY314/create_processed_tables.R"
dataset       <- "SDY314"
activity_name <- "process Immunespace, Import files into usable tables"
used          <- str_c(import_expr_id, import_trans_id, sep = ";") 

con   <- ImmuneSpaceR::CreateConnection(dataset)

ground_truth_df <-
    con$getDataset("fcs_analyzed_result") %>%
    dplyr::as_tibble() %>%
    mutate(sample = str_sub(participant_id, end = -5))

translation_df <- import_trans_id %>% 
    create_df_from_synapse_id() %>% 
    select("Subject Accession", "Expsample Accession") %>% 
    set_colnames(c("sample", "id"))

expr_df <- import_expr_id %>% 
    create_df_from_synapse_id() %>% 
    select(c("SYMBOL", contains("SIGNAL"))) %>% 
    dplyr::rename(Hugo = SYMBOL) %>% 
    gather(key = "id", value = "expr", -Hugo) %>% 
    drop_na() %>% 
    mutate(id = str_remove_all(id, ".AVG_SIGNAL")) %>% 
    inner_join(translation_df) %>% 
    select(sample, Hugo, expr)



# 
# ground_truth_immune_space <- 
#     con$getDataset("fcs_analyzed_result") %>% 
#     dplyr::as_tibble() %>% 
#     mutate(sample = str_sub(participant_id, end = -5))
# 
# ground_truth_10k <- whole_blood_flow_10k_id %>%
#     create_df_from_synapse_id %>%
#     filter(study_accession == dataset) %>%
#     dplyr::rename(sample = subject_accession) %>% 
#     select(sample, T_cells)

# samples_in_common <- intersect(expr_df$sample, ground_truth_immune_space$sample)
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
# # Is Illumina specific enough for the mirco array
# expression_manifest_df <- tibble(
#     path = c("expression_log.tsv", "expression_linear.tsv"),
#     parent = upload_id,
#     executed = script_url,
#     activityName = activity_name,
#     dataset = dataset,
#     used = used,
#     file_type = "expression",
#     expression_type = "microarray",
#     microarray_type = "Illumina",
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
