suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))
suppressPackageStartupMessages(p_load("data.table"))

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

exprs <-
    llply(synIds,
          .fun = function(synId) {
              obj <- synGet(synId, downloadFile = TRUE)
              read.table(obj$path, sep=",", header=TRUE, stringsAsFactors = FALSE)
          })

## These expression files have one sample that is not in the ground
## truth -- remove it.

synIds <-
    list("freq" = "syn22781220")

tbls <-
    llply(synIds,
          .fun = function(synId) {
              obj <- synGet(synId, downloadFile = TRUE)
              file <- obj$path
              read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)
          })

freq.mappings <-
    list("B.cells" = c("CD19..B.cells"),
        "CD4.T.cells" = c("CD4..T.cells"),
        "CD8.T.cells" = c("CD8..T.cells"),
        "NK.cells" = c("CD3..CD56."),
##        "neutrophils" = NA,
##        "monocytic.lineage" = NA,
##        "fibroblasts" = NA,
##        "endothelial.cells" = NA,
        "memory.B.cells" = c("Memory.B.cells"),
        "naive.B.cells" = c("Naive.B.cells"),
        "memory.CD4.T.cells" = c("CD4.central.memory", "CD4.EMRA", "CD4.EM"),
        "naive.CD4.T.cells" = c("CD4.naive"),
        "regulatory.T.cells" = c("CD4..T.reg"),
        "memory.CD8.T.cells" = c("CD8.central.memory", "CD8.EMRA", "CD8.EM"),
        "naive.CD8.T.cells" = c("CD8.naive"),
        "monocytes" = c("Monocytes")
##        "myeloid.dendritic.cells" = NA,
##        "macrophages" = NA
        )

translate.cell.types <- function(df, mapping, norm.col) {
    cols <- names(mapping)
    for(col in cols) {
        df[,col] <- unlist(apply(df[, mapping[[col]], drop=F], 1, function(row) sum(row)))
        df[,col] <- df[,col] / df[, norm.col]
    }
    df
}

extract.id <- function(str) {
    ret <- gsub(str, pattern=".fcs", replacement="")
    ret <- gsub(ret, pattern="Immune_", replacement="")
    ret <- gsub(ret, pattern=" Basal", replacement="")
    ret <- gsub(ret, pattern="Healthy BM ", replacement="")
    ret
}

freq.mappings <- lapply(freq.mappings, function(vec) paste0(vec, ".count"))
## norm.col <- "Live.cells.count"
## From Kara Davis:
## The library prep did not filter for live cells, that I know of. I would normalize by the “Cells.count” column as these are the cells that I gated as actual cells, not just debris.
norm.col <- "Cells.count"
print(colnames(tbls[["freq"]]))
translated.freqs <- translate.cell.types(tbls[["freq"]], freq.mappings, norm.col = norm.col)
translated.freqs$id <-
    unlist(lapply(translated.freqs$file, extract.id))

## These ids are numeric, which might create problems when they are the
## column names of the expression matrix. So append an "X", as will have
## been done implicitly by read.table for the expression matrices above.
translated.freqs$id <- paste0("X", translated.freqs$id)

expr.cols <- lapply(exprs, colnames)
expr.cols <- Reduce(intersect, expr.cols)

common.samples <- intersect(expr.cols, translated.freqs$id)

translated.freqs <- subset(translated.freqs, id %in% common.samples)

exprs <- llply(exprs, .fun = function(df) df[, c("Gene", common.samples)])

revised.files <-
    llply(files,
          .fun = function(file) gsub(file, pattern=".csv", replacement="-revised.csv"))

anno.tbl <- translated.freqs[, c("id", "Patient.ID", "Timepoint")]
flag <- is.na(anno.tbl$Timepoint)
anno.tbl[flag, "Timepoint"] <- "Normal"

translated.freqs <- translated.freqs[, c("id", names(freq.mappings))]
translated.freqs <- reshape2::melt(translated.freqs)
colnames(translated.freqs) <- c("sample.id", "cell.type", "measured")
translated.freqs <- cbind(dataset.name = "BALL", translated.freqs)

parent.id <- "syn22781213"

file <- "b-all-ground-truth.csv"
write.table(file = file, translated.freqs, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

file <- "b-all-annotation.csv"
write.table(file = file, anno.tbl, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

## Partition into coarse- and fine-grained 

## parent.id <- "syn22492020"
parent.id <- "syn23019061"

fine.cell.types <-
  c("memory.B.cells", "naive.B.cells", "memory.CD4.T.cells", "naive.CD4.T.cells", "regulatory.T.cells",
    "memory.CD8.T.cells", "naive.CD8.T.cells", "NK.cells", "neutrophils", "monocytes",
    "myeloid.dendritic.cells", "macrophages", "fibroblasts", "endothelial.cells")

coarse.cell.types <-
  c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage",
    "fibroblasts", "endothelial.cells")

fine.flag <- translated.freqs$cell.type %in% fine.cell.types
coarse.flag <- translated.freqs$cell.type %in% coarse.cell.types

if(!all(fine.flag | coarse.flag)) {
  print(translated.freqs[!(fine.flag | coarse.flag)])
  stop("Unknown cell type\n")
}

coarse.grained.gold.standard <- translated.freqs[coarse.flag, ]
fine.grained.gold.standard <- translated.freqs[fine.flag, ]
gold.standards <- list("ball-coarse-gold-standard.csv" =
                           coarse.grained.gold.standard,
                       "ball-fine-gold-standard.csv" =
                           fine.grained.gold.standard)


for(file in names(gold.standards)) {
    tbl <- gold.standards[[file]]
    write.table(file = file, tbl, sep=",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(paste0("Storing ", file, "\n"))
    f <- File(file, parentId = parent.id, synpaseStore = TRUE)
    synStore(f)
}


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
    "hugo.expr.file" = revised.files["hugo.expr.file"],
    "hugo.expr.est.counts.file" = revised.files["hugo.expr.est.counts.file"],
    "ensg.expr.file" = revised.files["ensg.expr.file"],
    "ensg.expr.est.counts.file" = revised.files["ensg.expr.est.counts.file"],
    "fastq.samples" = NA,                  
    "fastq1.files" = NA,
    "fastq2.files" = NA,
    "native.expr.file" = revised.files["native.expr.file"],
    "native.expr.est.counts.file" = revised.files["native.expr.est.counts.file"],
    "symbol.compression.est.counts.function" = "colMeans",
    "ensg.compression.est.counts.function" = "colMeans"
)
input.tbl <- as.data.frame(params)
fine.input.tbl <- input.tbl 
coarse.input.tbl <- fine.input.tbl

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

nms <- names(revised.files)
names(nms) <- nms
l_ply(nms,
      .fun = function(nm) {
          file <- revised.files[[nm]]
          expr <- exprs[[nm]]
          write.table(file = file, expr, row.names = FALSE, col.names = TRUE,
                      sep = ",", quote = FALSE)
          f <- File(file, parentId = parent.id, synapseStore = TRUE)
          synStore(f)          
      })

cat("Exiting successfully\n")
q(status=0)
