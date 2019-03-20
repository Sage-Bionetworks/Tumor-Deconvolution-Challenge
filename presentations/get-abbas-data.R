library(tidyverse)
library(synapser)
library(synapserutils)
library(data.table)
library(magrittr)
library(GEOquery)

gse.id <- "GSE22886"

gse_object <- getGEO(gse.id)

library(plyr)

nms <- names(gse_object)
names(nms) <- nms

source("../scripts/utils.R")

library(tools)

l_ply(nms,
  .fun = function(nm) {
           anno.tbl <- pData(phenoData(gse_object[[nm]]))
	   expr <- as.data.frame(exprs(gse_object[[nm]]))
	   gpl <- getGEO(gse_object[[nm]]@annotation, destdir=".")
	   mapping <- Table(gpl)[, c("ID", "Gene Symbol")]
	   colnames(mapping) <- c("from", "to")
	   if(!all(rownames(expr) %in% mapping$from)) {
	     cat("Some probes not in mapping\n")
	     print(table(rownames(expr) %in% mapping$from))
	     stop("Stopping")
	   } else {
             cat("All probes in mapping\n")
           }
	   expr_df <- expr %>%
	     compressGenes(mapping, fun = max) %>%
	     matrix_to_df("Hugo") %>%
             gather(key = "sample", value = "expr", - Hugo) %>%
             filter(Hugo != "") %>%
	     spread(key = "sample", value = "expr") %>%
	     as.data.frame()
           prefix <- file_path_sans_ext(nm)
	   write.table(file=paste0(prefix, "-anno.tsv"), anno.tbl, sep="\t", quote=FALSE,
	               row.names=FALSE, col.names=TRUE)
	   write.table(file=paste0(prefix, "-expr.tsv"), expr_df, sep="\t", quote=FALSE,
	               row.names=FALSE, col.names=TRUE)
		       
         })



