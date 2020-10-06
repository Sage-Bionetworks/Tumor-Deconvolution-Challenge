suppressPackageStartupMessages(library("pacman"))
suppressPackageStartupMessages(p_load("synapser"))
suppressPackageStartupMessages(p_load("plyr"))
suppressPackageStartupMessages(p_load("dplyr"))
suppressPackageStartupMessages(p_load("reshape2"))

synLogin()

synIds <-
    list("freq" = "syn22781220",
         "gated" = "syn22781219")

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
norm.col <- "Live.cells.count"
translated.freqs <- translate.cell.types(tbls[["freq"]], freq.mappings, norm.col = norm.col)
translated.freqs$id <-
    unlist(lapply(translated.freqs$file, extract.id))
translated.freqs <- translated.freqs[, c("id", names(freq.mappings))]
translated.freqs <- melt(translated.freqs)
colnames(translated.freqs) <- c("sample.id", "cell.type", "measured")
translated.freqs <- cbind(dataset.name = "BALL", translated.freqs)

parent.id <- "syn22781213"

file <- "b-all-ground-truth.csv"
write.table(file = file, translated.freqs, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, synapseStore = TRUE)
synStore(f)

## Partition into coarse- and fine-grained 

parent.id <- "syn22492020"

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

cat("Exiting successfully\n")
q(status=0)
