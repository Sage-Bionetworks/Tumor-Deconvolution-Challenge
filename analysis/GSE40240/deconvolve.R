library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

repo_dir  <- "../../../Tumor-Deconvolution-Challenge/"
tmp_dir   <- tempdir()
tmp_dir <- "."
## expr_id is GSEXXX/pre-processed/
## This expression file should be in linear/non-log space!
expr_id   <- "syn18162723"
## upload_id is GSEXXX/deconvolution_results/
upload_id <- "syn18163264"

dataset <- "GSE77343"
github.path <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
dataset.path <- paste0(github.path, dataset)
deconvolve.url <- paste0(dataset.path, "/deconvolve.R")

source(str_c(repo_dir, "scripts/synapseClient_functions.R"))
source(str_c(repo_dir, "scripts/matrix_functions.R"))
source(str_c(repo_dir, "scripts/deconvolution-utils.R"))
setwd(tmp_dir)
synapseLogin()

paths <- download.deconvolution.tools.and.matrices()
cs.cwltool <- paths$cs.cwltool
cs.ref.matrix <- paths$cs.ref.matrix
mcp.cwltool <- paths$mcp.cwltool

run.deconvolution.tools(expr_id, upload_id, cs.cwltool, cs.ref.matrix, mcp.cwltool, exec.url = deconvolve.url)
