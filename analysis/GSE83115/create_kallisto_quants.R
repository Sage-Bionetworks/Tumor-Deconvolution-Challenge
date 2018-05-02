library(tidyverse)
library(synapser)
library(doMC)

# local
# home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
# kallisto_dir <- "/home/aelamb/repos/kallisto_cwl/"

#ec2
home_dir  <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/ubuntu/tmp/"
kallisto_cwl_file <- "/home/ubuntu/kallisto_cwl/"


manifest_id          <- "syn12177468"
file_view_id         <- "syn12179146"
download_id_raw      <- "syn12177447"
download_id_sampled  <- "syn12177448"
upload_id            <- "syn12180291"


fasta_file <- "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
index_file <- "GRCH38.idx"    
script <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/SE83115/create_kallisto_quants.R"
cell_lines <- c("Monocytes_BGI", "CD4+_T_cells_BGI", "CD8+_T_cells_BGI")


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
n_cores = detectCores() - 1


manifest_df <- manifest_id %>% 
    create_df_from_synapse_id %>% 
    filter(sample_name %in% cell_lines)


file_view_df <- synTableQuery("select id, name, parentId from syn12179146") %>% 
    as.data.frame %>% 
    filter(str_detect(name, "fastq.gz$")) %>% 
    select(id, name, parentId)

raw_file_df <- file_view_df %>% 
    filter(parentId == download_id_raw) %>% 
    mutate(run = str_sub(name, end = 10)) %>% 
    inner_join(manifest_df) %>% 
    select(name, id, sample_name) %>% 
    mutate(pair = str_match(name, "^[:print:]+_([1,2]).fastq.gz$")[,2])


sampled_file_df <- file_view_df %>% 
    filter(parentId == download_id_sampled) %>% 
    mutate(pair = str_match(name, "^[:print:]+_p([1,2]).fastq.gz$")[,2]) %>% 
    mutate(sample_name = str_match(name, "^([:print:]+)_p[1,2].fastq.gz$")[,2])

file_df <- list(raw_file_df, sampled_file_df) %>% 
    bind_rows %>% 
    select(-parentId)
    


create_and_upload_quant_file <- function(df){
    fastq1 <- download_from_synapse(df$id[[1]])
    fastq2 <- download_from_synapse(df$id[[2]])
    cwl_command <- str_c("cwltool", 
                         str_c(kallisto_dir, "quant.cwl"),
                         "--index", index_file,
                         "--threads", n_cores,
                         "--fastq1", fastq1,
                         "--fastq2", fastq2, 
                         sep = " ")
    system(cwl_command)
    cwl_command2 <- str_c("cwltool", 
                          str_c(kallisto_dir, "h5dump.cwl"),
                          "--h5 abundance.h5",
                          sep = " ")
    system(cwl_command2)
    abundance_file <- str_c(df$sample_name[[1]], "_kallisto_abundance.tsv")
    file.rename("abundance.tsv", abundance_file)
    activity_obj <- Activity(
        name = "quantify and upload",
        description = "quantify fastq pair and upload to synapse",
        executed = list(script),
        used = flatten(list(manifest_id, file_view_id, fasta_file, df$id))
    )
    upload_file_to_synapse(abundance_file, upload_id, activity_obj = activity_obj)
    walk(c(abundance_file, "abundance.h5", "run_info.json"), file.remove)
}

file_df %>% 
    split(.$sample_name) %>% 
    walk(create_and_upload_quant_file)











