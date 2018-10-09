library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

repo_dir  <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir   <- tempdir()
expr_id   <- "syn13363401"
upload_id <- "syn13363404"


source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
cur_dir <- getwd()
setwd(tmp_dir)

synapseLogin()

print(tmp_dir)


expr_df <- expr_id %>% 
    create_df_from_synapse_id 

write_tsv(expr_df, "expr.tsv")

cs.cwltool <- "/home/aelamb/repos/irwg/iatlas-tool-cibersort/Dockstore.cwl"
cs.ref.matrix <- "/home/aelamb/repos/irwg/iatlas-tool-cibersort/sample.references.matrix.txt"
mcp.cwltool <- "/home/aelamb/repos/irwg/iatlas-tool-mcpcounter/Dockstore.cwl"

cs.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/Dockstore.cwl?token=AE4r2s5o0vfInYfTNeLoBE5dH68htLMyks5bxiiEwA%3D%3D"
cs.ref.matrix.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-cibersort/master/sample.references.matrix.txt?token=AE4r2nM3Ju0dq8GSfka8R-f-GLuCgpoCks5bxiiYwA%3D%3D"
mcp.cwltool.url <- "https://raw.githubusercontent.com/CRI-iAtlas/iatlas-tool-mcpcounter/master/Dockstore.cwl?token=AE4r2iAjO-QXRM5TTUR2vC7DXigoHUcGks5bxieHwA%3D%3D"

cs.cwltool <- "cs-Dockerstore.cwl"
cs.ref.matrix <- "sample.references.matrix.txt"
mcp.cwltool <- "mcp-Dockerstore.cwl"

urls <- c(cs.cwltool.url, cs.ref.matrix.url, mcp.cwltool.url)
files <- c(cs.cwltool, cs.ref.matrix, mcp.cwltool)
for(i in 1:length(urls)) {
    if(!file.exists(files[i])) {
        library(RCurl)
        download.file(urls[i], destfile=files[i], method="libcurl")
    }
}

download.file("https://d396qusza40orc.cloudfront.net/getdata%2Fdata%2Fss06hid.csv",destfile="reviews.csv",method="libcurl")

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
    paste("cwltool", mcp.cwltool),
    "--input_expression_file expr_matrix.tsv",
    "--output_file_string mcpcounter_results.tsv",
    "--features_type HUGO_symbols",
    sep = " "))

activity_obj <- Activity(
    name = "create",
    description = "create and upload deconvolution results using cibersort and mcpcounter cwl files",
    used = list(expr_id),
    executed = list("https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/GSE59654/fastq_mixing_CD4_CD8/deconvolve.R")
)

upload_file_to_synapse("cibersort_results.tsv", upload_id, activity_obj = activity_obj)
upload_file_to_synapse("mcpcounter_results.tsv", upload_id, activity_obj = activity_obj)

setwd(cur_dir)
unlink(tmp_dir, force = TRUE, recursive = TRUE)
