library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)


tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"
repo_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tpm_id    <- "syn16801064"
upload_id <- "syn16784234"


source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
setwd(tmp_dir)
synapseLogin()




tpm_df <- tpm_id %>%
    create_df_from_synapse_id()

write_tsv(tpm_df, "expr.tsv")

tpm_df %>% 
    df_to_matrix("Hugo") %>% 
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
    used = list(tpm_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE83115/fastq_mixing_CD4_CD8_2/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)
