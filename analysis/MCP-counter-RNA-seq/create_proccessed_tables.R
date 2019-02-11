library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(xlsx)

upload_id     <- "syn18207869"

fpkm_id         <- "syn18144910"
ground_truth_id <- "syn18144911"


source("../../scripts/utils.R")
synLogin()

dataset       <- "MCP-counter-RNA-seq"
github.path   <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path  <- paste0(github.path, dataset)
script_url    <- paste0(dataset.path, "/create_processed_tables.R")
used          <- str_c(fpkm_id, ground_truth_id, sep = ";")
activity_name <- "process data from MCPCounter data into usable tables"

fpkm_df <- fpkm_id %>% 
    create_df_from_synapse_id() %>% 
    dplyr::rename(Hugo = Gene.Symbol) %>% 
    gather(key = "sample", value = "expr", -"Hugo")

ground_truth_df <- ground_truth_id %>% 
    download_from_synapse() %>% 
    read.xlsx(sheetIndex = 1) %>% 
    dplyr::rename(cell_type = `NA.`) %>% 
    gather(key = "sample", value = "value", -"cell_type")
    

samples_in_common <- intersect(fpkm_df$sample, ground_truth_df$sample)


linear_expr_df <- fpkm_df %>%
    filter(sample %in% samples_in_common) %>%
    spread(key = "sample", value = "expr")

log_expr_df <- fpkm_df %>%
    filter(sample %in% samples_in_common) %>%
    mutate(expr = log10(expr + 1)) %>% 
    spread(key = "sample", value = "expr")

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common) %>%
    spread(key = "cell_type", value = "value")

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "expression",
    expression_type = "RNASeq",
    rnaseq_normalization = "FPKM",
    expression_space = c("log10", "linear")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    used = used,
    file_type = "ground truth",
    unit = "mcp score?",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")

syncToSynapse("expression_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")

