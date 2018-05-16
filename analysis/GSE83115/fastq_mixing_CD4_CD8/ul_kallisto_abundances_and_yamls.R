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

upload_id  <- "syn12231588"
yaml_file  <- "yaml.tsv"
fastq_file <- "fastq_files.tsv"

manifest_id  <- "syn12177468"
file_view_id <- "syn12179146"
index_id     <- "syn12213028"


setwd(home_dir)
source("scripts/utils.R")
setwd(work_dir)
synLogin()
registerDoMC(cores = detectCores() - 1)

yaml_df  <- yaml_file %>% 
    read_tsv %>% 
    mutate(tsv = str_c(str_remove(yaml, ".yaml"), ".tsv")) 

fastq_df <- read_tsv(fastq_file)

yaml_activity_obj <- Activity(
    name = "create and upload",
    description = "create yaml file for doing fastq mixing and kallisto quantification",
    executed = list(
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/dl_fastqs_for_mixing.R",
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/create_yamls_for_mixing_workflow.R",
        "https://github.com/Sage-Bionetworks/fastq_mixing_workflow_cwl"),
    used = c(manifest_id, file_view_id)
)
 
yaml_ids <- llply(
    yaml_df$yaml, 
    upload_file_to_synapse, 
    upload_id, 
    activity_obj = yaml_activity_obj, 
    ret = "syn_id", 
    .parallel = T)

yaml_df$yaml_id <- unlist(yaml_ids)

upload_tsvs_by_sample <- function(df){
    activity_obj <- Activity(
        name = "create and upload",
        description = "create yaml file for doing fastq mixing and kallisto quantification",
        executed = list(
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/dl_fastqs_for_mixing.R",
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/create_yamls_for_mixing_workflow.R",
            "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/ul_kallisto_abundances_and_yamls.R",
            "https://github.com/Sage-Bionetworks/fastq_mixing_workflow_cwl",
            "https://github.com/Sage-Bionetworks/fastq_mixer",
            "https://github.com/Sage-Bionetworks/kallisto_cwl"),
        used = c(df$yaml_id, fastq_df$id))
    id <- upload_file_to_synapse(df$tsv, upload_id, activity_obj = activity_obj, ret = "syn_id")
    return(id)
}


tsv_ids <- yaml_df %>% 
    split(1:nrow(.)) %>% 
    llply(upload_tsvs_by_sample, .parallel = F)

yaml_df$tsv_id <- unlist(tsv_ids)


yaml_df %>% 
    select(CD8_fractions, CD4_fractions, seed1, seed2, seed3, yaml, yaml_id, tsv, tsv_id) %>% 
    set_colnames(c("CD8_fraction", "CD4_fraction", "seed1", "seed2", "seed3", "yaml_name", "yaml_id", "abundance_name", "abundance_id")) %>% 
    write_tsv("manifest.tsv")

manifest_activity_obj <- Activity(
    name = "create and upload",
    description = "create yaml file for doing fastq mixing and kallisto quantification",
    executed = list(
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/dl_fastqs_for_mixing.R",
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/create_yamls_for_mixing_workflow.R",
        "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/fastq_mixing_CD4_CD8/ul_kallisto_abundances_and_yamls.R",
        "https://github.com/Sage-Bionetworks/fastq_mixing_workflow_cwl",
        "https://github.com/Sage-Bionetworks/fastq_mixer",
        "https://github.com/Sage-Bionetworks/kallisto_cwl"),
    used = c(yaml_df$yaml_id, fastq_df$id))    
    
upload_file_to_synapse("manifest.tsv", upload_id, activity_obj = manifest_activity_obj)










