library(tidyverse)
library(synapser)
library(doMC)

# local
# home_dir     <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# workflow_dir <- "/home/aelamb/repos/fastq_mixing_workflow_cwl/"
# work_dir     <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"

# ec2
home_dir     <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
workflow_dir <- "/home/ubuntu/fastq_mixing_workflow_cwl/"
work_dir     <- "/home/ubuntu/"

ground_truth_id     <- "syn12650217"
kallisto_index_id   <- "syn12213028"
synapse_config_file <- ".synapseConfig"
upload_id           <- "syn12678224"
p1_fastq_id         <- "syn12654356"
p2_fastq_id         <- "syn12656315"
donor_id            <- "HD30"
cell_type           <- "Monocytes"
min_readable_reads  <- 250000L
replicates          <- 5
multipliers         <- c(8, 4, 2, 1, 0.5, 0.25, 0.125)

setwd(home_dir)

source("scripts/utils.R")
setwd(work_dir)
synLogin()
n_cores <- detectCores() - 1


kallisto_threads <- as.integer(n_cores)

monocyte_p <- ground_truth_id %>% 
    create_df_from_synapse_id %>% 
    filter(`donor ID` == donor_id) %>% 
    extract2(cell_type)
    

# run specific parameters
yaml_df <-
    data_frame("multipliers" = multipliers) %>%
    mutate(mixer_total_reads = multipliers * min_readable_reads) %>% 
    mutate(mixer_total_reads = as.integer(mixer_total_reads)) %>% 
    merge(data_frame("rep" = 1:replicates)) %>% 
    mutate(prefix = str_c(as.character(mixer_total_reads), "_rep_", rep)) %>% 
    mutate(yaml = str_c(prefix, ".yaml")) %>% 
    inset("mixer_seed", value = sample(1:10000, nrow(.)))

source(str_c(workflow_dir, "utils/write_yaml.R"))

# create yamls
create_synapse_workflow_yaml_by_row <- function(row){
    create_synapse_workflow_yaml(
        yaml_file = row$yaml,
        synapse_config_file = synapse_config_file,
        output_name = row$prefix,
        fastq_p1_synapse_ids = p1_fastq_id,
        fastq_p2_synapse_ids = p2_fastq_id,
        mixer_fractions = 1.0,
        upload_id = upload_id,
        kallisto_index_synapse_id = kallisto_index_id,
        mixer_seed = row$mixer_seed,
        mixer_total_reads = row$mixer_total_reads,
        kallisto_threads = kallisto_threads,
        annotations = list(
            "seed" = row$mixer_seed,
            "run" = row$rep,
            "reads" = row$mixer_total_reads))
}

yaml_df %>% 
    split(1:nrow(.)) %>% 
    walk(create_synapse_workflow_yaml_by_row)

