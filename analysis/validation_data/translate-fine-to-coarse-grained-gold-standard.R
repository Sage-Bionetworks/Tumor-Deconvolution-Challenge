suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(synapser))

source("../utils.R")

synLogin()

synId <- "syn21820376"

file <- synGet(synId, downloadFile=TRUE)$path
fine.grained.gold.standard.tbl <- read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)

translate.fine.grained.gold.standard <- function(fine.grained.gold.standard.tbl) {

    tmp <- acast(fine.grained.gold.standard.tbl , sample.id ~ cell.type, value.var = "measured", fill = NA)

    ## B.cells = memory.B.cells + naive.B.cells
    flag <- colnames(tmp) %in% c("memory.B.cells", "naive.B.cells")
    tmp <- cbind(tmp, B.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row, na.rm=TRUE)))))

    ## CD4.T.cells = memory.CD4.T.cells + naive.CD4.T.cells + regulatory.T.cells
    flag <- colnames(tmp) %in% c("memory.CD4.T.cells", "naive.CD4.T.cells", "regulatory.T.cells")
    tmp <- cbind(tmp, CD4.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row, na.rm=TRUE)))))
    
    ## CD8.T.cells = memory.CD8.T.cells + naive.CD8.T.cells
    flag <- colnames(tmp) %in% c("memory.CD8.T.cells", "naive.CD8.T.cells")
    tmp <- cbind(tmp, CD8.T.cells = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row, na.rm=TRUE)))))
    
    ## monocytic.lineage = 
    flag <- colnames(tmp) %in% c("myeloid.dendritic.cells", "macrophages", "monocytes")
    tmp <- cbind(tmp, monocytic.lineage = as.numeric(unlist(apply(tmp[, flag, drop=F], 1, function(row) sum(row, na.rm=TRUE)))))
    
    tmp <- tmp[, colnames(tmp) %in% coarse.cell.types]
    
    coarse.grained.gold.standard <- melt(tmp)
    colnames(coarse.grained.gold.standard) <- c("sample.id", "cell.type", "measured")

    coarse.grained.gold.standard

}

trans.tbl <-
    ddply(fine.grained.gold.standard.tbl,
          .variables = c("dataset.name"),
          .fun = function(df) {
              translate.fine.grained.gold.standard(df)
          })

## Upload this translated ground truth to synapse
folder.synId <- "syn21820011"
file <- "coarse-translated-from-fine.csv"
write.table(file = file, trans.tbl, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

obj <- File(file, parent = folder.synId)
synStore(obj)

## Also, create the input file (by subsetting the one used for the challenge to just the fine-grained datasets)
synId <- "syn21821122"
file <- synGet(synId, downloadFile=TRUE)$path
input.tbl <- read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)

revised.input.tbl <- subset(input.tbl, dataset.name %in% trans.tbl$dataset.name)

folder.synId <- "syn21821096"
file <- "input-translated-from-fine.csv"
write.table(file = file, revised.input.tbl, sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)

obj <- File(file, parent = folder.synId)
synStore(obj)
