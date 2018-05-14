library(plyr)
library(doMC)
library(tidyverse)
library(magrittr)
library(synapser)

# local
# home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# work_dir <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"


#ec2
home_dir <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
work_dir <- "/home/ubuntu/"

setwd(home_dir)
source("scripts/utils.R")
setwd(work_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

upload_id  <- "syn12180291"
fastq_file <- "fastq2.tsv"

manifest_id          <- "syn12177468"
file_view_id         <- "syn12179146"

manifest_df <- read_tsv(fastq_file) 

yaml_activity_obj <- Activity(
    name = "create and upload",
    description = "create yaml file for doing kallisto quantification",
    executed = list(
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/dl_fastqs_for_kallisto.R",
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/create_yamls_for_kallisto_quant.R",
        "https://github.com/Sage-Bionetworks/kallisto_cwl",
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/ul_kallisto_abundances_and_yamls.R"),
    used = c(manifest_id, file_view_id)
)
 
ids <- llply(
    manifest_df$yaml, 
    upload_file_to_synapse, 
    upload_id, 
    activity_obj = yaml_activity_obj, 
    ret == "syn_id", 
    .parallel = T)

manifest_df$yaml_id <- ids

upload_tsvs_by_sample <- function(df){
    activity_obj <- Activity(
        name = "create and upload",
        description = "kallisto quantification",
        executed = list(
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/dl_fastqs_for_sampling.R",
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/create_yamls_for_sampling.R",
            "https://github.com/Sage-Bionetworks/fastq_mixer",
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_sampling_analysis/ul_sampled_fastqs_and_yamls.R"),
        used = c(manifest_id, file_view_id, df$p1_id, df$p2_id, df$yaml_id))
    id <- upload_file_to_synapse(df$tsv, upload_id, activity_obj = activity_obj, return == "syn_id")
    return(id)
}

ids <- manifest_df %>% 
    split(.$yaml) %>% 
    llply(upload_tsvs_by_sample, .parallel = F)

manifest_df$tsv_id <- ids

manifest_df %>% 
    select(sample_name, p1_id, p2_id, yaml, yaml_id, tsv, tsv_id) %>% 
    set_colnames(c("sample_name", "fastq_p1_id", "fastq_p2_id", "yaml_name", "yaml_id", "abundance_name", "asundance_id")) %>% 
    write_tsv("manifest.tsv")

upload_file_to_synapse("manifest.tsv", upload_id, activity_obj = yaml_activity_obj)










