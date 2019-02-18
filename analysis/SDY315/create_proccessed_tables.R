library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(Biobase)
library(ImmuneSpaceR)

upload_id     <- "syn18347648"

source("../../scripts/utils.R")
synLogin()

script_url    <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY315/create_processed_tables.R"
dataset       <- "SDY315"
activity_name <- "process Immunespace files into usable tables"

con   <- ImmuneSpaceR::CreateConnection(dataset)
expr_set <- con$getGEMatrix(c("SDY315_Other_GroupA_Geo", "SDY315_Other_GroupB_Geo", "SDY315_Other_GroupC_Geo"))

expr_df <- con$mapSampleNames(EM = expr_set) %>% 
    Biobase::exprs() %>% 
    matrix_to_df("Hugo") %>% 
    gather(key = "sample", value = "expr",  -Hugo) %>% 
    mutate(sample = str_remove(sample, ".315_d0"))

ground_truth_df <- 
    con$getDataset("fcs_analyzed_result") %>% 
    dplyr::as_tibble() %>% 
    mutate(sample = str_sub(participant_id, end = -5))


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
# ground_truth_df <- ground_truth_df %>%
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
