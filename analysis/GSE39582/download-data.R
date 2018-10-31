suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(library("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

library(GEOquery)
gse <- getGEO("GSE39582", GSEMatrix=TRUE)

ex <- as.data.frame(exprs(gse[[1]]))

file <- "GSE39582-expr.tsv"
write.table(file = file, ex, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

gpl <- getGEO(gse[[1]]@annotation, destdir=".")
mapping <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(mapping) <- c("from", "to")
if(!all(rownames(ex) %in% mapping$from)) {
    cat("Some probes not in mapping\n")
    table(rownames(ex) %in% mapping$from)
    stop("Stopping")
} else {
    cat("All probes in mapping\n")
}

library(plyr)
## Translate/compress genes from one name space (e.g., probe ids) to another (e.g., symbols)
compressGenes <- function(e, mapping, from.col = "from", to.col = "to")
{
  e$to    <- mapping[match(rownames(e), mapping[, from.col]), to.col]
  e           <- e[!is.na(e$to),]
  e           <- ddply(.data = e, .variables = "to", .fun = function(x){apply(x[,-ncol(x)],2,mean)},.parallel = T)
  rownames(e) <- e$to
  e           <- e[,-1]
  return(e)
}

ex.translated <- compressGenes(ex, mapping)

file <- "GSE39582-expr-symbol.tsv"
write.table(file = file, ex.translated, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
