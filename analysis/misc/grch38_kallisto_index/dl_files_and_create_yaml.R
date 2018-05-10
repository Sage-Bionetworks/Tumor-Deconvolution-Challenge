library(tidyverse)
library(magrittr)

# local
# home_dir     <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
# kallisto_dir <- "/home/aelamb/repos/kallisto_cwl/"
# work_dir     <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

# ec2
home_dir     <- "/home/ubuntu/Tumor-Deconvolution-Challenge/"
kallisto_dir <- "/home/ubuntu/kallisto_cwl/"
work_dir     <- "/home/ubuntu/"

fasta_url     <- "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
kallisto_yaml <- "grch38_kallisto_index.yaml"
index_file    <- "GRCH38.idx"

setwd(home_dir)
source("scripts/utils.R")
source(str_c(kallisto_dir, "utils/write_yaml.R"))
setwd(work_dir)
synLogin()

system(str_c("wget ", fasta_url))
fasta_file <- basename(fasta_url)

create_kallisto_index_yaml(kallisto_yaml, fasta_file, index_file)
