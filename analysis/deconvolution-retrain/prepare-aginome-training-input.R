suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(optparse))
suppressPackageStartupMessages(p_load(data.table))

option_list <- list(
     make_option(c("--input-csv-file"), type="character", default=NULL,
                 help="Input CSV expression matrix of Challenge data", metavar="character"),
     make_option(c("--output-purified-tsv-file"), type="character", default=NULL,
                 help="Output TSV expression matrix of _purified_ Challenge data", metavar="character"),
#     make_option(c("--output-admixtures-tsv-file"), type="character", default=NULL,
#                 help="Output TSV expression matrix of _admixtures_ from Challenge data", metavar="character"),
     make_option(c("--output-purified-csv-anno-file"), type="character", default=NULL,
                 help="Output CSV annotation file of _purified_ Challenge data", metavar="character"),
     make_option(c("--grain"), type="character", default=NULL,
                 help="coarse or fine", metavar="character")

                 )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

source("../utils.R")

input.csv.file <- opt$`input-csv-file`
output.tsv.file <- opt$`output-purified-tsv-file`
#output.admixtures.tsv.file <- opt$`output-admixtures-tsv-file`
output.csv.anno.file <- opt$`output-purified-csv-anno-file`

input.mat <- as.data.frame(fread(input.csv.file))
rownames(input.mat) <- input.mat$Gene

trans.tbl <- get.purified.sample.translation.table()
rownames(trans.tbl) <- trans.tbl$sample

purified.flag <- (colnames(input.mat) %in% rownames(trans.tbl))
out.mat <- input.mat[, purified.flag]
#admixture.flag <- !(colnames(input.mat) %in% rownames(trans.tbl)) & !(colnames(input.mat) %in% c("CRC", "Breast", "Gene"))
#out.admixtures.mat <- input.mat[, admixture.flag]

trans.col <- NULL
if(opt$grain == "coarse") {
  trans.col <- "coarse.grained.cell.type"
} else if(opt$grain == "fine") {
  trans.col <- "fine.grained.cell.type"
}

if(is.null(trans.col)) { stop("stop") }

common <- intersect(colnames(input.mat), rownames(trans.tbl))
input.mat <- input.mat[, common]
trans.tbl <- trans.tbl[common, ]

trans.tbl <- trans.tbl[, c("sample", trans.col)]
colnames(trans.tbl) <- c("sample.name", "cell.type")

write.table(file = output.csv.anno.file, trans.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")

write.table(file = output.tsv.file, out.mat, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

output.csv.file <- gsub(output.tsv.file, pattern=".tsv", replacement=".csv")
write.table(file = output.csv.file, out.mat, row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")

#write.table(file = output.admixtures.tsv.file, out.admixtures.mat, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
#output.admixtures.csv.file <- gsub(output.admixtures.tsv.file, pattern=".tsv", replacement=".csv")
#write.table(file = output.admixtures.csv.file, out.admixtures.mat, row.names=TRUE, col.names=TRUE, quote=FALSE, sep=",")