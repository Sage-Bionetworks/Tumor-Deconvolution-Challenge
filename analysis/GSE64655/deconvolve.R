library(plyr)
library(doMC)
library(tidyverse)
library(synapser)
library(data.table)
library(magrittr)

home_dir <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir  <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"

expr_id      <- "syn11969378"
hugo_id       <- "syn11536071"

upload_id  <- "syn11969392"


setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()
registerDoMC(cores = detectCores())

hugo_df <-  create_df_from_synapse_id(hugo_id)

expr_df <- expr_id %>%
    create_df_from_synapse_id(unzip = T, skip = 3) %>% 
    select(-c(`Gene Type`, Description, `Gene Symbol`)) %>% 
    left_join(hugo_df, by = c("Gene ID" = "ensembl_gene_id")) %>% 
    select(-`Gene ID`) %>% 
    group_by(hgnc_symbol) %>% 
    summarise_all(sum) %>% 
    ungroup %>% 
    filter(hgnc_symbol != "")

write_tsv(expr_df, "expr.tsv")

expr_df %>% 
    df_to_matrix("hgnc_symbol") %>% 
    write.table("expr_matrix.tsv", sep = "\t", quote = F)

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-cibersort/Dockstore.cwl", 
    "--mixture_file expr.tsv", 
    "--sig_matrix_file /home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt", 
    "--QN",
    "--output_file_string cibersort_results.tsv",
    sep = " "))

system(str_c(
    "cwltool /home/aelamb/repos/irwg/iatlas-tool-mcpcounter/Dockstore.cwl",
    "--input_expression_file expr_matrix.tsv",
    "--output_file_string mcpcounter_results.tsv",
    "--features_type HUGO_symbols",
    sep = " "))

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(expr_id, hugo_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE64655/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)



