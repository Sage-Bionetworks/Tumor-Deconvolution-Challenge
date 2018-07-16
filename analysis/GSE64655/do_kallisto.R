library(doMC)
library(plyr)
library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)


# home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"
# cwl_file <- "/home/aelamb/repos/kallisto_cwl/fastq_abundances_workflow.cwl"

home_dir <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/ubuntu/tmp/"
cwl_file <- "/home/ubuntu/kallisto_cwl/fastq_abundances_workflow.cwl"

manifest_id <- "syn12663606"
file_dir_id <- "syn12649849"
upload_id   <- "syn13841771" 
index_id    <- "syn12213028"

synapseCacheDir(tmp_dir)
registerDoMC(cores = 4)
synapseClient::synapseLogin()
setwd(tmp_dir)

source(str_c(home_dir, "scripts/synapseClient_functions.R"))

index_file <- download_from_synapse(index_id)

file_df <- file_dir_id %>% 
    get_file_df_from_synapse_dir_id %>% 
    filter(str_detect(file.name, ".fastq.gz$")) %>% 
    mutate(SRR_id = str_match(file.name, "(SRR[0-9]*)_[0-9]{1}.fastq.gz")[,2]) %>% 
    mutate(pair = str_match(file.name, "SRR[0-9]*_([0-9]{1}).fastq.gz")[,2]) %>% 
    select(-file.name) %>% 
    spread(key = "pair", value = file.id) %>% 
    set_colnames(c("SRR_id", "synapse_id1", "synapse_id2"))

manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(day == 0) %>% 
    left_join(file_df) %>% 
    select(cell_type, patient, synapse_id1, synapse_id2) %>% 
    mutate(cell_type = str_replace_all(cell_type, " ", "_")) %>% 
    mutate(output_name = str_c(patient, cell_type, "abundance.tsv", sep = "_")) 

write_tsv(manifest_df, "manifest.tsv")

activity_obj <- Activity(
    name = "run kallisto",
    used = list(manifest_id, file_dir_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/do_kallisto.R"))

upload_file_to_synapse("manifest.tsv", upload_id, activity_obj = activity_obj)

do_kallisto_by_row <- function(args){
    fastq_file1 <- download_from_synapse(args$synapse_id1)
    fastq_file2 <- download_from_synapse(args$synapse_id2)
    command <- str_c(
        "cwltool", 
        "cwl_file",
        "--index_file", index_file,
        "--fastq_file1", fastq_file1,
        "--fastq_file2", fastq_file2)
    system(command)
    activity_obj <- Activity(
        name = "run kallisto",
        used = list(args$synapse_id1, args$synapse_id2, index_id),
        executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/do_kallisto.R"))
    file.rename("out/abundance.tsv", args$output_name)
    upload_file_to_synapse(args$output_name, upload_id, activity_obj = activity_obj)
    file.remove(fastq_file1)
    file.remove(fastq_file2)
}

manifest_df %>% 
    split(1:nrow(.)) %>% 
    l_ply(do_kallisto_by_row, .parallel = T)





