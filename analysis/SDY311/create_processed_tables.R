library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(tidyverse)

upload_id     <- "syn18347671"

fc_wb_10K_id      <- "syn13363370"
cytof_pbmc_10K_id <- "syn13363372"
import_cytof_id   <- "syn18347685"
import_expr_id    <- "syn18347682"
import_trans_id   <- "syn18347684"


source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY311/create_processed_tables.R"
dataset       <- "SDY311"
activity_name <- "process Import 10K files into usable tables"
used          <- str_c(fc_wb_10K_id, cytof_pbmc_10K_id, import_cytof_id, import_expr_id, import_trans_id, sep = ";")



ground_truth_flow_wb_10K_df <- fc_wb_10K_id %>%
    create_df_from_synapse_id %>%
    filter(study_accession == dataset) %>%
    .[, colSums(is.na(.)) < nrow(.)] %>%
    dplyr::rename(sample = subject_accession)

ground_truth_cytof_pbmc_10K_df <- cytof_pbmc_10K_id %>%
    create_df_from_synapse_id %>%
    filter(study_accession == dataset) %>% 
    dplyr::rename(sample = subject_accession)

ground_truth_cytof_import <- import_cytof_id %>% 
    create_df_from_synapse_id %>% 
    dplyr::rename(sample = `Subject Accession`)

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
    mutate(id = str_remove_all(id, "_SIGNAL")) %>% 
    inner_join(translation_df) %>% 
    select(sample, Hugo, expr)


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
