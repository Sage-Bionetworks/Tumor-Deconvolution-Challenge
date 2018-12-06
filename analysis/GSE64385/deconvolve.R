library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

repo_dir  <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir   <- tempdir()
tmp_dir <- "."
## expr_id is GSEXXX/pre-processed/<expr-file>.tsv
expr_id   <- "syn15667753"
## upload_id is GSEXXX/deconvolution_results/
upload_id <- "syn17088610"


source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
setwd(tmp_dir)
synapseLogin()

cs.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/Dockstore.cwl?token=AE4r2kVqI149gsoBR0izIReUJrMPUVwkks5b14pZwA%3D%3D"
cs.ref.matrix.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/sample.references.matrix.txt?token=AE4r2jd69jitwyGzi7oMfmYSrvkRqraoks5b14xPwA%3D%3D"
mcp.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-mcpcounter/master/Dockstore.cwl?token=AE4r2iOJ9rWrdgPc5ogjS1k850vvK6I8ks5b14ynwA%3D%3D"

cs.cwltool <- "../../external/cs-Dockerstore.cwl"
cs.ref.matrix <- "../../external/sample.references.matrix.txt"
mcp.cwltool <- "../../external/mcp-Dockerstore.cwl"

urls <- c(cs.cwltool.url, cs.ref.matrix.url, mcp.cwltool.url)
files <- c(cs.cwltool, cs.ref.matrix, mcp.cwltool)
for(i in 1:length(urls)) {
  if(!file.exists(files[i])) {
    library(RCurl)
    cat(paste0("Downloading ", files[i], "\n"))
    download.file(urls[i], destfile=files[i], method="libcurl")
  }
}

expr_df <- expr_id %>% 
    create_df_from_synapse_id 

write_tsv(expr_df, "expr.tsv")

expr_df %>% 
    df_to_matrix("Hugo") %>% 
    write.table("expr_matrix.tsv", sep = "\t", quote = F)

system(str_c(
    paste("cwltool", cs.cwltool),
    "--mixture_file expr.tsv", 
    paste("--sig_matrix_file", cs.ref.matrix),
    "--output_file_string cibersort_results.tsv",
    sep = " "))

system(str_c(
    paste("cwltool", cs.cwltool),
    "--mixture_file expr.tsv", 
    paste("--sig_matrix_file", cs.ref.matrix),
    "--output_file_string cibersort_abs_results.tsv",
    "--abs_method no.sumto1",
    "--absolute",
    sep = " "))

system(str_c(
    paste("cwltool", mcp.cwltool),
    "--input_expression_file expr_matrix.tsv",
    "--output_file_string mcpcounter_results.tsv",
    "--features_type HUGO_symbols",
    sep = " "))

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE65133/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("cibersort_abs_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)
