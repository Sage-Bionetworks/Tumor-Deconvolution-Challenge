library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

source("../../scripts/utils.R")
synLogin()

expression_id <- "syn17014418"
ground_truth_id <- "syn17014449"
upload_id <- "syn17022196"

tempdir <- "/home/aelamb/tmp/tumor_deconvolution/dry_run_files/"

gt_df_long <- ground_truth_id %>% 
    create_df_from_synapse_id() %>% 
    gather(key = "cell_type", value = "value", - sample) 

expr_df_long <- expression_id %>% 
    create_df_from_synapse_id() %>% 
    gather(key = "sample", value = "value", - Hugo) 

samples_in_common <- intersect(gt_df_long$sample, expr_df_long$sample)

gt_df <- gt_df_long %>% 
    filter(sample %in% samples_in_common) %>% 
    mutate(cell_type = str_sub(cell_type, end = -5)) %>% 
    spread(key = "sample", value = "value")

expr_df <- expr_df_long %>% 
    filter(sample %in% samples_in_common) %>% 
    separate(Hugo, into = "Gene", sep = " /// ", extra = "drop") %>% 
    group_by(Gene, sample) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    spread(key = "sample", value = "value")

remove(gt_df_long, expr_df_long)

write_csv(gt_df, str_c(tempdir, "goldstandard.csv"))

activity_obj1 <- Activity(
    name = "create file",
    description = "goldstandard for dryrun",
    used = list(ground_truth_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/challenge/dry_run/create_dry_run_files.R"))

upload_file_to_synapse(str_c(tempdir, "goldstandard.csv"), upload_id, activity_obj = activity_obj1)


write_csv(expr_df, str_c(tempdir, "input.csv"))

activity_obj2 <- Activity(
    name = "create file",
    description = "goldstandard for dryrun",
    used = list(expression_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/challenge/dry_run/create_dry_run_files.R"))

upload_file_to_synapse(str_c(tempdir, "input.csv"), upload_id, activity_obj = activity_obj2)
    