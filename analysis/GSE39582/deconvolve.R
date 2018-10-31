library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

repo_dir  <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir   <- tempdir()
tmp_dir <- "."

## leaderboard_datasets/GSE39572/pre-processed/expression_microarray.tsv
expr_id   <- "syn17014418"

## leaderboard_datasets/GSE39572/deconvolution_results/
upload_id <- "syn17014457"

source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
cur_dir <- getwd()
setwd(tmp_dir)

synapseLogin()

expr_df <- expr_id %>% 
    create_df_from_synapse_id 

write_tsv(expr_df, "expr.tsv")


cs.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/Dockstore.cwl?token=AE4r2kVqI149gsoBR0izIReUJrMPUVwkks5b14pZwA%3D%3D"
cs.ref.matrix.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/sample.references.matrix.txt?token=AE4r2jd69jitwyGzi7oMfmYSrvkRqraoks5b14xPwA%3D%3D"
mcp.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-mcpcounter/master/Dockstore.cwl?token=AE4r2iOJ9rWrdgPc5ogjS1k850vvK6I8ks5b14ynwA%3D%3D"

cs.cwltool <- "cs-Dockerstore.cwl"
cs.ref.matrix <- "sample.references.matrix.txt"
mcp.cwltool <- "mcp-Dockerstore.cwl"

urls <- c(cs.cwltool.url, cs.ref.matrix.url, mcp.cwltool.url)
files <- c(cs.cwltool, cs.ref.matrix, mcp.cwltool)
for(i in 1:length(urls)) {
    if(!file.exists(files[i])) {
        library(RCurl)
        cat(paste0("Downloading ", files[i], "\n"))
        download.file(urls[i], destfile=files[i], method="libcurl")
    }
}

expr_df %>% 
    df_to_matrix("Hugo") %>% 
    write.table("expr_matrix.tsv", sep = "\t", quote = F)

source("CIBERSORT.R")

cat(str_c(
    paste("cwltool", cs.cwltool),
    "--mixture_file expr.tsv", 
    paste("--sig_matrix_file", cs.ref.matrix),
    "--output_file_string cibersort_results.tsv",
    sep = " "))

cs.matrix <- read.table(cs.ref.matrix, sep="\t", header=TRUE, row.names=1)

stop("stop")

## result_matrix <- CIBERSORT(cs.ref.matrix, "expr.tsv", perm=0, QN=TRUE, absolute=TRUE, abs_method='no.sumto1')

result_matrix <- CIBERSORT_(cs.matrix, expr, perm=0, QN=TRUE, absolute=TRUE, abs_method='no.sumto1')

ofile <- "cibersort_results.tsv"

result_matrix %>% 
    as.data.frame %>% 
    rownames_to_column("sample") %>% 
    write_tsv(ofile)

stop("stop")


system(str_c(
    paste("cwltool", cs.cwltool),
    "--mixture_file expr.tsv", 
    paste("--sig_matrix_file", cs.ref.matrix),
    "--output_file_string cibersort_results.tsv",
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
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE39582/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)

setwd(cur_dir)
if(tmp_dir != ".") { unlink(tmp_dir, force = TRUE, recursive = TRUE) }
