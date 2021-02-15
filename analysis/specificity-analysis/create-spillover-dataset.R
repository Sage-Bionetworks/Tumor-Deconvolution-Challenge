suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))

synLogin()

## Read in the expression data (of different formats) and subset
## to the purified samples.
purified.samples <- c(
    "Naive_B_cells_1",
    "Macrophages_2",
    "Dendritic_cells_1",
    "CRC",
    "Macrophages_1",
    "Fibroblasts",
    "Endothelial_cells",
    "Memory_CD8_T_cells_1",
    "Monocytes_1",
    "Neutrophils_2",
    "Breast",
    "Dendritic_cells_2",
    "NK_cells_2",
    "Monocytes_2",
    "NK_cells_1",
    "Memory_CD8_T_cells_2",
    "Memory_CD4_T_cells_2",
    "Naive_CD4_T_cells_1",
    "Tregs",
    "Memory_CD4_T_cells_1",
    "Naive_CD8_T_cells_2",
    "Naive_CD4_T_cells_2")

cols.to.keep <- c("Gene", purified.samples)

samples <- purified.samples
names(samples) <- samples

## Map populations to samples
populations <- samples
populations <- gsub(populations, pattern="_1", replacement="")
populations <- gsub(populations, pattern="_2", replacement="")
pop.df <- data.frame(population = unname(populations), sample = names(populations))

## Create the gold standard files for the purified samples
## Translate the population names to those we use in the challenge
fine.grained.map <- list(
    "Breast" = "Breast",
    "CRC" = "CRC",
    "Dendritic_cells" = "myeloid.dendritic.cells",
    "Endothelial_cells" = "endothelial.cells",
    "Fibroblasts" = "fibroblasts",
    "Macrophages" = "macrophages",
    "Memory_CD4_T_cells" = "memory.CD4.T.cells",
    "Memory_CD8_T_cells" = "memory.CD8.T.cells",
    "Monocytes" = "monocytes",
    "NK_cells" = "NK.cells",
    "Naive_B_cells" = "naive.B.cells",
    "Naive_CD4_T_cells" = "naive.CD4.T.cells",
    "Naive_CD8_T_cells" = "naive.CD8.T.cells",
    "Neutrophils" = "neutrophils",
    "Tregs" = "regulatory.T.cells")

fine.grained.map.df <- data.frame(population = names(fine.grained.map), challenge.population = unname(unlist(fine.grained.map)))
fine.grained.pop.df <- merge(pop.df, fine.grained.map.df)

## Add purified samples to the gold standard
tmp <- fine.grained.pop.df
tmp$val <- 1
tmp <- tmp[, c("sample", "challenge.population", "val")]
tmp2 <- melt(acast(tmp, sample ~ challenge.population, fill = 0))
colnames(tmp2) <- c("sample.id", "cell.type", "measured")
tmp2 <- subset(tmp2, (cell.type != "Breast") & (cell.type != "CRC"))
gold.standard.tbl <- cbind(dataset.name = "DS5", tmp2)
for(col in c("dataset.name", "sample.id", "cell.type")) {
    gold.standard.tbl[, col] <- as.character(gold.standard.tbl[, col])
}

## fine-grained
# memory.B.cells (missing)
# naive.B.cells
# memory.CD4.T.cells
# naive.CD4.T.cells
# regulatory.T.cells
# memory.CD8.T.cells
# naive.CD8.T.cells
# NK.cells
# neutrophils (missing)
# monocytes
# myeloid.dendritic.cells
# macrophages
# fibroblasts
# endothelial.cells
fine.grained.gold.standard <- gold.standard.tbl

## coarse-grained
# B.cells
# CD4.T.cells
# CD8.T.cells
# NK.cells
# neutrophils (missing)
# monocytic.lineage
# fibroblasts
# endothelial.cells
coarse.cell.types <-
  c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage",
    "fibroblasts", "endothelial.cells")

tmp <- acast(gold.standard.tbl , sample.id ~ cell.type, value.var = "measured", fill = 0)

## B.cells = memory.B.cells + naive.B.cells
flag <- colnames(tmp) %in% c("memory.B.cells", "naive.B.cells")
tmp <- cbind(tmp, B.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## CD4.T.cells = memory.CD4.T.cells + naive.CD4.T.cells + regulatory.T.cells
flag <- colnames(tmp) %in% c("memory.CD4.T.cells", "naive.CD4.T.cells", "regulatory.T.cells")
tmp <- cbind(tmp, CD4.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## CD8.T.cells = memory.CD8.T.cells + naive.CD8.T.cells
flag <- colnames(tmp) %in% c("memory.CD8.T.cells", "naive.CD8.T.cells")
tmp <- cbind(tmp, CD8.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

## monocytic.lineage = 
flag <- colnames(tmp) %in% c("myeloid.dendritic.cells", "macrophages", "monocytes")
tmp <- cbind(tmp, monocytic.lineage = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row)))))

tmp <- tmp[, colnames(tmp) %in% coarse.cell.types]

coarse.grained.gold.standard <- melt(tmp)

colnames(coarse.grained.gold.standard) <- c("sample.id", "cell.type", "measured")

coarse.grained.gold.standard <- cbind(dataset.name = "DS5",
                                      coarse.grained.gold.standard)

gold.standards <- list("specificity-coarse-gold-standard.csv" =
                           coarse.grained.gold.standard,
                       "specificity-fine-gold-standard.csv" =
                           fine.grained.gold.standard)


purified.matrices <- list()

## this is symbol_tpm.csv
synId <- "syn21576632"
obj <- synGet(synId, downloadFile = TRUE)
cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

purified.matrices[["symbol_tpm"]] <-
    list("file" = "symbol_tpm.csv", "mat" = cpm.expr[, cols.to.keep])

## this is ensg_tpm.csv
synId <- "syn21576631"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cpm.expr <- read.table(obj$path, sep = ",", header = TRUE)

purified.matrices[["ensg_tpm"]] <-
    list("file" = "ensg_tpm.csv", "mat" = ensg.cpm.expr[, cols.to.keep])

## this is symbol_counts.csv
synId <- "syn21576630"
obj <- synGet(synId, downloadFile = TRUE)
cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

purified.matrices[["symbol_counts"]] <-
    list("file" = "symbol_counts.csv", "mat" = cnts.expr[, cols.to.keep])

## this is ensg_counts.csv
synId <- "syn21576629"
obj <- synGet(synId, downloadFile = TRUE)
ensg.cnts.expr <- read.table(obj$path, sep = ",", header = TRUE)

purified.matrices[["ensg_counts"]] <-
    list("file" = "ensg_counts.csv", "mat" = ensg.cnts.expr[, cols.to.keep])

ds <- "DS5"
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
    "hugo.expr.file" = purified.matrices[["symbol_tpm"]][["file"]],
    "hugo.expr.est.counts.file" = purified.matrices[["symbol_counts"]][["file"]],
    "ensg.expr.file" = purified.matrices[["ensg_tpm"]][["file"]],
    "ensg.expr.est.counts.file" = purified.matrices[["ensg_counts"]][["file"]],
    "fastq.samples" = NA,                  
    "fastq1.files" = NA,
    "fastq2.files" = NA,
    "native.expr.file" = purified.matrices[["ensg_tpm"]][["file"]],
    "native.expr.est.counts.file" = purified.matrices[["ensg_counts"]][["file"]],
    "symbol.compression.est.counts.function" = "colMeans",
    "ensg.compression.est.counts.function" = "colMeans"
)
fine.input.tbl <- as.data.frame(params)
coarse.input.tbl <- fine.input.tbl

parent.id <- "syn22392130"

for(nm in names(purified.matrices)) {
    file <- purified.matrices[[nm]][["file"]]
    mat <- purified.matrices[[nm]][["mat"]]
    write.table(file = file, mat, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
    cat(paste0("Storing ", file, "\n"))
    f <- File(file, parentId = parent.id, synpaseStore = TRUE)
    synStore(f)
    file.remove(file)
}
    
for(file in names(gold.standards)) {
    tbl <- gold.standards[[file]]
    write.table(file = file, tbl, sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(paste0("Storing ", file, "\n"))
    f <- File(file, parentId = parent.id, synpaseStore = TRUE)
    synStore(f)
}

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
