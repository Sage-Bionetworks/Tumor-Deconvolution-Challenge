library(plyr)
library(doMC)
library(tidyverse)
library(magrittr)
library(synapser)

# local
home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
work_dir <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"


#ec2
# home_dir <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
# work_dir <- "/home/ubuntu/"

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
    mutate(prefix = str_remove(yaml, ".yaml")) %>% 
    mutate(output_name1 = str_c(prefix, "_rep1.tsv")) %>%
    mutate(output_name2 = str_c(prefix, "_rep2.tsv")) %>% 
    mutate(output_name3 = str_c(prefix, "_rep3.tsv")) %>% 
    select(- c(log, prefix))
              
fastq_df <- fastq_file %>% 
    read_tsv %>% 
    select(-path)

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
    id1 <- upload_file_to_synapse(df$output_name1, upload_id, activity_obj = activity_obj, ret = "syn_id")
    id2 <- upload_file_to_synapse(df$output_name2, upload_id, activity_obj = activity_obj, ret = "syn_id")
    id3 <- upload_file_to_synapse(df$output_name3, upload_id, activity_obj = activity_obj, ret = "syn_id")
    return(c(id1, id2, id3))
}


output_ids <- yaml_df %>% 
    split(1:nrow(.)) %>% 
    llply(upload_tsvs_by_sample, .parallel = T)

yaml_df$output_id1 <- map_chr(output_ids, function(ids) ids[[1]])
yaml_df$output_id2 <- map_chr(output_ids, function(ids) ids[[2]])
yaml_df$output_id3 <- map_chr(output_ids, function(ids) ids[[3]])


yaml_df %>% 
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










