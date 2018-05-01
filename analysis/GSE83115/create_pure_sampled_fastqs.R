library(tidyverse)
library(synapser)

# local
# home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

#ec2
home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/ubuntu/tmp/"

manifest_id  <- "syn12177468"
file_view_id <- "syn12179146"
download_id  <- "syn12177447"
upload_id    <- "syn12177448"

script <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/create_pure_sampled_fastqs.R"

cell_lines_to_sample <- c("Monocytes_BGI", 
                          "CD4+_T_cells_BGI",
                          "CD8+_T_cells_BGI")

setwd(home_dir)
source("scripts/utils.R")
source("scripts/combine_fastqs.R")
setwd(tmp_dir)
synLogin()

manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines_to_sample)

file_view_df <- synTableQuery("select id, name, parentId from syn12179146") %>% 
    as.data.frame %>% 
    filter(parentId == download_id) %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    mutate(run = str_sub(name, end = 10)) %>% 
    inner_join(manifest_df) %>% 
    select(name, id, sample_name)

create_pure_samples <- function(df){
    paths <- map(df$id, download_from_synapse)
    paths %>%
        str_c("gunzip ", .) %>%
        walk(system)
    new_paths <- str_sub(paths, end = - 4)
    combine_df <- data_frame(
        "fraction" = c(1),
        "p1_fastq_file" = c(new_paths[[1]]),
        "p2_fastq_file" = c(new_paths[[2]]))
    output_files <- combine_paired_fastq_files_by_n(combine_df, 3, df$sample_name[[1]])
    walk(new_paths, file.remove)
    output_files %>% 
        str_c("gzip ", .) %>% 
        walk(system)
    gzip_output_files <- str_c(output_files, ".gz")
    activity_obj <- Activity(
        name = "sample and upload",
        description = "sample with replacement from original fastq",
        executed = list(script),
        used = flatten(list(manifest_id, file_view_id, df$id))
    )
    walk(gzip_output_files, upload_file_to_synapse, upload_id, activity_obj)
    walk(gzip_output_files, file.remove)
}

df <- file_view_df %>% 
    split(.$sample_name) %>% 
    walk(create_pure_samples)
    


