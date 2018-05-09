library(plyr)
library(doMC)
library(tidyverse)
library(magrittr)
library(synapser)

# local
# home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# work_dir <- "/home/aelamb/"


#ec2
home_dir <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
work_dir <- "/home/ubuntu/"

setwd(home_dir)
source("scripts/utils.R")
setwd(work_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

upload_id    <- "syn12177448"
manifest_id  <- "syn12177468"
file_view_id <- "syn12179146"
download_id  <- "syn12177447"

cell_lines_to_sample <- c("Monocytes_BGI", 
                          "CD4+_T_cells_BGI",
                          "CD8+_T_cells_BGI")



manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    select(run, sample_name) %>% 
    set_names(c("SRR_id", "cell_line")) %>% 
    filter(cell_line %in% cell_lines_to_sample) %>%
    mutate(cell_line = str_remove(cell_line, "\\+"))

file_view_df <- file_view_id %>% 
    str_c("select id, name, parentId from ", .) %>% 
    synTableQuery %>% 
    as.data.frame %>% 
    filter(parentId == download_id) %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    mutate(SRR_id = str_sub(name, end = 10)) %>% 
    mutate(sample = str_sub(name, end = 12)) %>% 
    mutate(pair = str_sub(sample, start = 12)) %>% 
    select(id, SRR_id, pair)  %>% 
    spread(key = pair, value = id) %>% 
    set_names(c("SRR_id", "p1_source_id", "p2_source_id")) %>% 
    inner_join(manifest_df)

file_df <- list.files() %>% 
    keep(str_detect(., ".yaml$")) %>% 
    data_frame("yaml" = .) %>% 
    mutate(sample = str_match(yaml, "([:print:]+).yaml$")[,2]) %>% 
    mutate(cell_line = str_match(sample, "([:print:]+)_output[0-9]{1}$")[,2]) %>% 
    mutate(p1_fastq = str_c(sample, "_p1.fastq")) %>% 
    mutate(p2_fastq = str_c(sample, "_p2.fastq")) %>% 
    mutate(p1_fastq_gz = str_c(p1_fastq, ".gz")) %>% 
    mutate(p2_fastq_gz = str_c(p2_fastq, ".gz")) %>% 
    inner_join(file_view_df)



# 
# yaml_activity_obj <- Activity(
#     name = "create and upload",
#     description = "create yaml file for doing fastq sampling",
#     executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/ul_sampled_fastqs_and_yamls.R",
#                     "https://github.com/Sage-Bionetworks/fastq_mixer"),
#     used = c(manifest_id, file_view_id)
# )
# 
# 
# l_ply(file_df$yaml, upload_file_to_synapse, upload_id, activity_obj = yaml_activity_obj, .parallel = T)
# 
# gzip_commands <- str_c("gzip ", c(file_df$p1_fasta, file_df$p2_fasta))
# l_ply(gzip_commands, system, .parallel = T) 
# 
# upload_fastqs_by_sample <- function(df){
#     
#     
#     activity_obj <- Activity(
#         name = "create and upload",
#         description = "sample from fastq files",
#         executed = list(
#             "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/dl_fastqs_for_sampling.R",
#             "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/create_yamls_for_sampling.R",
#             "https://github.com/Sage-Bionetworks/fastq_mixer",
#             "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/ul_sampled_fastqs_and_yamls.R"),
#         used = c(manifest_id, file_view_id))
# }










