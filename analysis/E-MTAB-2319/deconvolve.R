library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)
library(MCPcounter)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/E-MTAB-2319/"

count_id      <- "syn11958709"
hugo_id       <- "syn11536071"

upload_id  <- "syn11968722"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

hugo_df <-  create_df_from_synapse_id(hugo_id)

cpm_df <- count_id %>%
    create_df_from_synapse_id(unzip = T) %>% 
    select(-c(Chr, Start, End, Strand, Length)) %>% 
    df_to_matrix("Geneid") %>% 
    .[rowSums(.) > 0,] %>%  
    add(1) %>% 
    apply(2, calculate_cpm) %>% 
    matrix_to_df("ensembl_gene_id") %>% 
    inner_join(hugo_df) %>% 
    select(hgnc_symbol, everything()) %>% 
    filter(hgnc_symbol != "") %>% 
    select(-ensembl_gene_id) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(.funs = sum) %>% 
    ungroup

write_tsv(cpm_df, "cpm.tsv")

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-cibersort/Dockstore.cwl", 
    "--mixture_file cpm.tsv", 
    "--sig_matrix_file /home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt", 
    "--QN",
    "--output_file_string cibersort_results.tsv",
    sep = " "))

cpm_df %>% 
    df_to_matrix("hgnc_symbol") %>% 
    write.table("cpm_matrix.tsv", sep = "\t", quote = F)

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-mcpcounter/Dockstore.cwl",
    "--input_expression_file cpm_matrix.tsv",
    "--output_file_string mcpcounter_results.tsv",
    "--features_type HUGO_symbols",
    sep = " "))

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(count_id, hugo_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/E-MTAB-2319/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)



