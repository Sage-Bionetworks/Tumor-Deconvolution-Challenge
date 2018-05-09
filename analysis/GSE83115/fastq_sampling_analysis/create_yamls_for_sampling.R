library(plyr)
library(doMC)
library(tidyverse)
library(yaml)

# local
# home_dir         <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir          <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
# cache_dir        <- "/home/aelamb/.synapseCache/"
# fastq_mixer_repo <- "/home/aelamb/repos/fastq_mixer/"

#ec2
home_dir         <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir          <- "/home/ubuntu/tmp/"
cache_dir        <- "/home/ubuntu/.synapseCache/"
fastq_mixer_repo <- "/home/ubuntu/fastq_mixer/"
   
manifest_id  <- "syn12177468"
file_view_id <- "syn12179146"
download_id  <- "syn12177447"


cell_lines_to_sample <- c("Monocytes_BGI", 
                          "CD4+_T_cells_BGI",
                          "CD8+_T_cells_BGI")

setwd(home_dir)
source("scripts/utils.R")
setwd(fastq_mixer_repo)
source("utils/write_yaml.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines_to_sample) %>% 
    select(run, sample_name) %>% 
    rename("SRR_id" = run)

fastq_df <- cache_dir %>% 
    list.files(recursive = T) %>% 
    keep(str_detect(., ".fastq$")) %>% 
    data_frame("path" = .) %>% 
    mutate(file_name = basename(path)) %>% 
    mutate(SRR_id = str_match(file_name, "([:alnum:]+)_")[,2]) %>% 
    inner_join(manifest_df) %>% 
    mutate(pair = str_match(file_name, "[:alnum:]+_([12]{1}).fastq$")[,2]) %>% 
    arrange(sample_name, pair)

create_yaml_by_sample <- function(df){
    map(1:3, function(run_number) create_yaml(df, run_number))
}

create_yaml <- function(df, seed){
    output_prefix  <- str_c(fastq_df$sample_name[[1]], "_output", as.character(seed))
    
    create_fastq_mixer_yaml(
        yaml_file = str_c( output_prefix , ".yaml"),
        fastq_files_p1 = fastq_df$path[[1]], 
        fastq_files_p2 = fastq_df$path[[2]], 
        sample_fractions = 1.0,
        seed = seed,
        output_prefix = output_prefix)
}

fastq_df %>% 
    split(.$sample_name) %>% 
    walk(create_yaml_by_sample)

