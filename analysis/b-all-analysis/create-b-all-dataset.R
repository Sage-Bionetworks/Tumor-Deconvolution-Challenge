suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))

synLogin()

hugo.tpm.synId <- "syn22492068"
hugo.est.counts.synId <- "syn22492053"
ensg.tpm.synId <- "syn22492055"
ensg.est.counts.synId <- "syn22492047"

synIds <-
  list("hugo.expr.file" = hugo.tpm.synId,
       "hugo.expr.est.counts.file" = hugo.est.counts.synId,
       "ensg.expr.file" = ensg.tpm.synId,
       "ensg.expr.est.counts.file" = ensg.est.counts.synId,
       "native.expr.file" = ensg.tpm.synId,
       "native.expr.est.counts.file" = ensg.est.counts.synId)

files <-
    llply(synIds,
          .fun = function(synId) {
              obj <- synGet(synId, downloadFile = FALSE)
              obj$properties$name
          })

ds <- "BALL"
ct <- NA
params <- list(
    "dataset.name" = ds,
    "cancer.type" = ct,
    "platform" = "Illumina",
    "scale" = "Linear",
    "normalization" = "TPM",
    "native.probe.type" = "ENSG",
    "symbol.compression.function" = "colMeans",
    "ensg.compression.function" = "colMeans",
    "symbol.to.native.mapping.file" = "native_to_hugo.tsv",
    "ensg.to.native.mapping.file" = "native_to_ensg.tsv",
    "hugo.expr.file" = files["hugo.expr.file"],
    "hugo.expr.est.counts.file" = files["hugo.expr.est.counts.file"],
    "ensg.expr.file" = files["ensg.expr.file"],
    "ensg.expr.est.counts.file" = files["ensg.expr.est.counts.file"],
    "fastq.samples" = NA,                  
    "fastq1.files" = NA,
    "fastq2.files" = NA,
    "native.expr.file" = files["native.expr.file"],
    "native.expr.est.counts.file" = files["native.expr.est.counts.file"],
    "symbol.compression.est.counts.function" = "colMeans",
    "ensg.compression.est.counts.function" = "colMeans"
)
input.tbl <- as.data.frame(params)
fine.input.tbl <- input.tbl 
coarse.input.tbl <- fine.input.tbl

parent.id <- "syn22492020"

file <- "input.csv"
write.table(file = file, input.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "fine-input.csv"
write.table(file = file, fine.input.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "coarse-input.csv"
write.table(file = file, coarse.input.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

cat("Exiting successfully\n")
q(status=0)
