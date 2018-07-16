library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"
cwl_file <- "/home/aelamb/repos/kallisto_cwl/fastq_abundances_workflow.cwl"

manifest_id <- "syn12663606"
file_dir_id <- "syn12649849"
upload_id   <- "syn13841771" 

synapseClient::synapseLogin()
setwd(tmp_dir)

source(str_c(home_dir, "scripts/synapseClient_functions.R"))

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
    select(cell_type, patient, synapse_id1, synapse_id2)

upload_file_to_synapse("out/abundance.tsv", upload_id, activity_obj = activity_obj)

do_kallisto_by_row <- function(args){
    command <- str_c(
        "cwltool", 
        "cwl_file",
        "--index_file", index_file,
        "--threads", n_cores,
        "--fastq_file1", args$synapse_id1,
        "--fastq_file2", args$synapse_id2)
    system(command)
    activity_obj <- Activity(
        name = "run kallisto",
        used = list(args$synapse_id1, args$synapse_id2, index_id),
        executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/do_kallisto.R"))
    upload_file_to_synapse("out/abundance.tsv", upload_id, activity_obj = activity_obj)

}





