library(synapser)


home_dir  <- "/home/aelamb/repos/Tumor-Deconvolution-Challenge/"
tmp_dir   <- "/home/aelamb/tmp/tumor_deconvolution/GSE83115/"

upload_id  <- "syn12212845"
index_file <- "GRCH38.idx"

# used
fasta_url     <- "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
kallisto_yaml <- "grch38_kallisto_index.yaml"

# executed
script1 <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/misc/grch38_kallisto_index/dl_files_and_create_yaml.R"
script2 <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/misc/grch38_kallisto_index/ul_index_and_yaml.R"
cwl     <- "https://github.com/Sage-Bionetworks/kallisto_cwl/kallisto_index.cwl"

setwd(home_dir)
source("scripts/utils.R")
setwd(tmp_dir)
synLogin()


yaml_activity_obj <- Activity(
    name = "create and upload",
    description = "create GRCH38 kallisto index",
    executed = list(script1, script2, cwl),
    used = c(fasta_url)
)

yaml_id <- upload_file_to_synapse(
    kallisto_yaml,
    upload_id,
    activity_obj = yaml_activity_obj,
    return = "syn_id"
)

index_activity_obj <- Activity(
    name = "create and upload",
    description = "create GRCH38 kallisto index",
    executed = list(script1, script2, cwl),
    used = c(fasta_url, yaml_id)
)

upload_file_to_synapse(
    index_file,
    upload_id,
    activity_obj = index_activity_obj
)


