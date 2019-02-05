library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(ImmuneSpaceR)
library(Biobase)
library(xlsx)

gt_excel_sheet <- "syn18070332"
upload_id     <- "syn18070255"
gt_upload_id  <- "syn18070253"

source("../../scripts/utils.R")
synLogin()

dataset <- "SDY67"
script_url <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SDY67/create_processed_tables.R"
activity_name <- "Process files from ImmuneSpacee."
used <- gt_excel_sheet


study <- ImmuneSpaceR::CreateConnection(dataset, verbose=TRUE)
expressionset1 <- study$getGEMatrix("SDY67_PBMC_HealthyAdults")
expressionset2 <- study$mapSampleNames(EM = expressionset1)

sample_names <- expressionset2 %>% 
    pData() %>% 
    filter(study_time_collected == 0) %>% 
    mutate(sample = str_sub(participant_id, end = -4)) %>% 
    use_series(sample) %>% 
    unique %>% 
    sort

expr_df <- expressionset2 %>% 
    exprs %>% 
    matrix_to_df("Hugo") %>% 
    arrange(Hugo) %>% 
    gather(key = "sample", value = "expr", - Hugo) %>%
    separate(col = sample, sep = "[[:punct:]]", into = c("sample", "x", "day")) %>% 
    filter(day == "d0") %>% 
    select(-c(x, day)) %>% 
    filter(sample %in% sample_names)




    

ground_truth_df <- gt_excel_sheet %>% 
    download_from_synapse() %>% 
    read.xlsx(2, startRow = 4) %>% 
    dplyr::rename(sample = donor.ID) %>% 
    gather(key = "cell_type", value = "percent", - sample) %>% 
    mutate(cell_type = str_replace_all(cell_type, "\\.", "_")) %>% 
    mutate(fraction = percent / 100) %>% 
    select(-percent)

samples_in_common <- intersect(ground_truth_df$sample, expr_df$sample)

tpm_expr_m <- expr_df %>%
    filter(sample %in% samples_in_common) %>% 
    mutate(expr = 2^expr) %>% 
    spread(key = "sample", value = "expr") %>% 
    df_to_matrix("Hugo") %>% 
    calculate_cpm() 

log_expr_df <- tpm_expr_m %>%
    add(1) %>% 
    log10() %>% 
    matrix_to_df("Hugo")

linear_expr_df <- matrix_to_df(tpm_expr_m, "Hugo")

ground_truth_df <- ground_truth_df %>%
    filter(sample %in% samples_in_common) %>% 
    spread(key = "cell_type", value = "fraction")
    

write_tsv(linear_expr_df, "expression_linear.tsv")
write_tsv(log_expr_df, "expression_log.tsv")
write_tsv(ground_truth_df, "ground_truth.tsv")

expression_manifest_df <- tibble(
    path = c("expression_log.tsv", "expression_linear.tsv"),
    parent = upload_id,
    used = used,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "expression",
    expression_type = "RNASeq", 
    rnaseq_normalization = "TPM",
    expression_space = c("log10", "linear")
)

ground_truth_manifest_df <- tibble(
    path = "ground_truth.tsv",
    parent = gt_upload_id,
    executed = script_url,
    activityName = activity_name,
    dataset = dataset,
    file_type = "ground truth",
    unit = "fraction",
    cell_types = str_c(colnames(ground_truth_df)[-1], collapse = ";")
)

write_tsv(expression_manifest_df, "expression_manifest.tsv")
write_tsv(ground_truth_manifest_df, "ground_truth_manifest.tsv")


syncToSynapse("expression_manifest.tsv")
syncToSynapse("ground_truth_manifest.tsv")
