library(synapseClient)
library(doMC)


upload_id       <- "syn12678224"
p1_fastq_id     <- "syn12654356"
p2_fastq_id     <- "syn12656315"
# kallisto_index_id   <- "syn12213028"
cwl_index_path  <- "/home/aelamb/repos/kallisto_cwl/index.cwl"
cwl_wf_path     <- "/home/aelamb/repos/kallisto_cwl/fastq_abundances_workflow.cwl"
work_dir        <- "/home/aelamb/tmp/tumor_deconvolution/GSE64655/"
fasta_url       <- "ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"


synapseLogin()

setwd(work_dir)



p1_fastq_path  <- synGet(p1_fastq_id)@filePath
p2_fastq_path  <- synGet(p2_fastq_id)@filePath
threads        <- detectCores() - 1

system(paste("wget", fasta_url))


cwl_command1 <- paste(
    "cwltool",
    cwl_index_path,
    "--index_file_string", "Grch38.idx",
    "--fasta_file", basename(fasta_url))

system(cwl_command1)


cwl_command2 <- paste(
    "cwltool",
    cwl_wf_path,
    "--index_file", "Grch38.idx",
    "--threads", threads,
    "--fastq_file1", p1_fastq_path,
    "--fastq_file2", p2_fastq_path)

system(cwl_command2)

entity <- File(path = "abundance.tsv", parentId = upload_id, synapseStore = TRUE)
synStore(entity)

