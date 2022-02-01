
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(xlsx))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(ggcorrplot))
suppressPackageStartupMessages(p_load(cowplot))
suppressPackageStartupMessages(p_load(edgeR)) # for cpm

synLogin()

use.log <- FALSE

figs.dir <- "figs/"
dir.create(figs.dir)

# Evidently CD45RO is a post-transcriptional modification? Don't show this transcript.
include.cd45ro <- FALSE
if(include.cd45ro) {
  ## Get the _transcript_ expression matrix
  synId <- "syn21574261"
  trans.expr.mat.file <- synGet(synId, downloadFile = TRUE)$path
  trans.expr.mat <- read.table(trans.expr.mat.file, sep="\t", header=TRUE)

  ## CD45RO = ENST00000348564.11 = PTPRC-201
  cd45ro.enst.id <- "ENST00000348564.11"
  rownames(trans.expr.mat) <- as.character(trans.expr.mat$transcript)
  trans.expr.mat <- trans.expr.mat[, !(colnames(trans.expr.mat) == "transcript")]
}

## Get the pseudo-count expression matrix
synId <- "syn21576630"
expr.mat.file <- synGet(synId, downloadFile = TRUE)$path

cat("Reading expression matrix\n")
expr.mat <- as.data.frame(fread(expr.mat.file))
rownames(expr.mat) <- as.character(expr.mat$Gene)
expr.mat <- expr.mat[, !(colnames(expr.mat) == "Gene")]
cat("Done reading expression matrix\n")

## Only keep the purified samples
colnames(expr.mat)[colnames(expr.mat) == "Breast"] <- "BRCA"
if(include.cd45ro) {
  colnames(trans.expr.mat)[colnames(trans.expr.mat) == "Breast"] <- "BRCA"
}
cols <- colnames(expr.mat)
##cols <- cols[!(grepl(cols, pattern="BM"))]
##cols <- cols[!(grepl(cols, pattern="RM"))]

sample.levels <- c(
    "Naive_B_cells_1",
    "Memory_CD4_T_cells_1",
    "Memory_CD4_T_cells_2",        
    "Naive_CD4_T_cells_1",
    "Naive_CD4_T_cells_2",    
    "Tregs",
    "Memory_CD8_T_cells_1",
    "Memory_CD8_T_cells_2",        
    "Naive_CD8_T_cells_2",
    "NK_cells_1",
    "NK_cells_2",        
    "Neutrophils_2",
    "Monocytes_1",
    "Monocytes_2",        
    "Dendritic_cells_1",
    "Dendritic_cells_2",           
    "Macrophages_1",
    "Macrophages_2",        
    "Endothelial_cells",
    "Fibroblasts",
    "BRCA",
    "CRC"
)
cols <- sample.levels

expr.mat <- expr.mat[, cols]
if(include.cd45ro) {
  trans.expr.mat <- trans.expr.mat[, cols]
}

if(include.cd45ro) {
  expr.mat <- rbind(expr.mat, CD45RO = as.numeric(trans.expr.mat[cd45ro.enst.id,]))
}

prioritize.top.genes <- function(log.expr.mat, n.top = 1000) {
  # make sure this first cutoff is done in log space, even
  # if we otherwise use linear space
  cutoff <- log2(1)
  if(!use.log) { cutoff <- 2 }
  keep.genes <- unlist(apply(log.expr.mat, 1, function(row) any(row > cutoff)))
  means <- unlist(apply(log.expr.mat, 1, mean))
  sds <- unlist(apply(log.expr.mat, 1, sd))
  cvs <- sds / means
  cvs <- cvs[keep.genes]
  cvs <- cvs[order(cvs, decreasing = TRUE)]
  keep.genes <- names(cvs)[1:n.top]
  return(keep.genes)
}

challenge.expr.mat <- expr.mat

challenge.expr.mat <- 
  plyr::rename(challenge.expr.mat, 
    replace = c("Naive_B_cells_1" = "Naive B",
                "Memory_CD4_T_cells_1" = "Memory CD4 T 1",
		"Memory_CD4_T_cells_2" = "Memory CD4 T 2",
		"Naive_CD4_T_cells_1" = "Naive CD4 T 1",
		"Naive_CD4_T_cells_2" = "Naive CD4 T 2",
		"Memory_CD8_T_cells_1" = "Memory CD8 T 1",
		"Memory_CD8_T_cells_2" = "Memory CD8 T 2",
		"Naive_CD8_T_cells_2" = "Naive CD8 T",
		"NK_cells_1" = "NK 1",
		"NK_cells_2" = "NK 2",
		"Neutrophils_2" = "Neutrophils",
		"Monocytes_1" = "Mono 1",
		"Monocytes_2" = "Mono 2",
		"Dendritic_cells_1" = "DCs 1",
		"Dendritic_cells_2" = "DCs 2",
		"Macrophages_1" = "Macro 1",
		"Macrophages_2" = "Macro 2",
		"Endothelial_cells" = "Endo",
		"Fibroblasts" = "Fibro"))




# Read in the Abbas data

abbas.expr.file <- paste0("../../presentations/","GSE22886-GPL96_series_matrix.txt-expr.tsv")
abbas.anno.file <- paste0("../../presentations/","GSE22886-GPL96_series_matrix.txt-anno.tsv")

process.abbas.data <- function(abbas.expr.file, abbas.anno.file) {

  orig.expr <- read.table(abbas.expr.file, sep="\t", header=TRUE)
  orig.anno <- read.table(abbas.anno.file, sep="\t", header=TRUE, as.is=TRUE)

  expr <- orig.expr
  anno <- orig.anno

  common.samples <- intersect(orig.anno$geo_accession, colnames(orig.expr))
  anno <- subset(orig.anno, geo_accession %in% common.samples)
  common.samples <- anno$geo_accession

  expr <- orig.expr[, c("Hugo", common.samples)]

  cell.type.means <- 
      dlply(anno, .variables = c("cell.type.ch1"),
            .fun = function(anno.cell) {
                expr.sub <- expr[, colnames(expr) %in% anno.cell$geo_accession]
                tmp <- rowMeans(expr.sub)
                if(use.log) { tmp <- log2(tmp) }
                tmp
            })
  tmp <- ldply(cell.type.means)
  rownames(tmp) <- tmp$cell.type.ch1
  tmp <- tmp[, !(colnames(tmp) == "cell.type.ch1")] 
  cell.type.means <- t(tmp)

  old.cols <- c("B cells", "IgM memory B cells", "CD8+ T cells",
                "CD4+ T cells", "CD4+ CD45RO+ CD45RA- T cells",
                "NK cells")
  new.cols <- c("B cells", "Memory B cells", "CD8 T cells",
                "CD4 T cells", "Memory CD4 cells", "NK cells")
  cell.type.means <- cell.type.means[, old.cols]
  colnames(cell.type.means) <- new.cols

  mat <- cell.type.means[, c("B cells", "CD8 T cells", "CD4 T cells", "Memory CD4 cells", "NK cells")]
  rownames(mat) <- expr$Hugo
  return(mat)
}




cibersort.supp.xls.file <- "41592_2015_BFnmeth3337_MOESM207_ESM.xls"

process.cibersort.data <- function(cibersort.supp.xls.file) {
  tmp <- read.xlsx(cibersort.supp.xls.file, sheetIndex=1, startRow=14)
  rownames(tmp) <- tmp$Gene.symbol
  tmp <- tmp[, 4:ncol(tmp)]
  if(use.log) {
    tmp <- log2(tmp)
  }
  return(tmp)
}

cibersort.log.expr.mat <- process.cibersort.data(cibersort.supp.xls.file)

n.top <- nrow(cibersort.log.expr.mat)
cibersort.genes <- rownames(cibersort.log.expr.mat)
n.top <- 1000
cat(paste0("n.top = ", n.top, "\n"))

replace = c("B.cells.naive" = "Naive B",
	    "B.cells.memory" = "Memory B",
            "T.cells.CD8" = "CD8 T",
            "T.cells.CD4.naive" = "Naive CD4 T",
            "T.cells.CD4.memory.resting" = "Memory Resting CD4 T",
            "T.cells.CD4.memory.activated" = "Memory Activated CD4 T",
            "T.cells.regulatory..Tregs." = "Tregs",
            "NK.cells.resting" = "Resting NK",
            "NK.cells.activated" = "Activated NK",
            "Monocytes" = "Mono",
            "Macrophages.M0" = "Macro M0",
            "Macrophages.M1" = "Macro M1",
            "Macrophages.M2" = "Macro M2",
            "Dendritic.cells.resting" = "Resting DCs",
            "Dendritic.cells.activated" = "Activated DCs",
            "Neutrophils" = "Neutrophils")

cibersort.log.expr.mat <-
  plyr::rename(cibersort.log.expr.mat, replace = replace)

cibersort.log.expr.mat <- cibersort.log.expr.mat[, unname(replace)]

cibersort.cor.mat <- cor(cibersort.log.expr.mat, method="pearson")

challenge.log.expr.mat <- cpm(challenge.expr.mat, log=use.log)
challenge.top.genes <- prioritize.top.genes(challenge.log.expr.mat, n.top = n.top)
#challenge.top.genes <- rownames(challenge.log.expr.mat)
#challenge.top.genes <- challenge.top.genes[challenge.top.genes %in% cibersort.genes]
challenge.log.expr.mat <- challenge.log.expr.mat[challenge.top.genes,]
challenge.cor.mat <- cor(challenge.log.expr.mat, method = "pearson")

abbas.log.expr.mat <- process.abbas.data(abbas.expr.file, abbas.anno.file)
colnames(abbas.log.expr.mat) <- gsub(colnames(abbas.log.expr.mat), pattern=" cells", replacement="")
abbas.top.genes <- prioritize.top.genes(abbas.log.expr.mat, n.top = n.top)
#abbas.top.genes <- rownames(abbas.log.expr.mat)
#abbas.top.genes <- abbas.top.genes[abbas.top.genes %in% cibersort.genes]
abbas.log.expr.mat <- abbas.log.expr.mat[abbas.top.genes,]
abbas.cor.mat <- cor(abbas.log.expr.mat, method = "pearson")

cat(paste0("abbas num genes: ", nrow(abbas.log.expr.mat), "\n"))
cat(paste0("challenge num genes: ", nrow(challenge.log.expr.mat), "\n"))
cat(paste0("cibersort num genes: ", nrow(cibersort.log.expr.mat), "\n"))

g.challenge <- ggcorrplot(challenge.cor.mat, hc.order = TRUE, hc.method = "average")
g.abbas <- ggcorrplot(abbas.cor.mat, hc.order = TRUE, hc.method = "average")
g.cibersort <- ggcorrplot(cibersort.cor.mat, hc.order = TRUE, hc.method = "average")

file.prefix <- paste0(figs.dir, "/", "signature-correlations")
png(paste0(file.prefix, ".png"), width = 2 * 480, height = 2 * 480)
g.all <- plot_grid(g.challenge, g.abbas, g.cibersort, nrow=2, ncol=2, labels="AUTO")
print(g.all)
d <- dev.off()

write.table(challenge.cor.mat, file = paste0(figs.dir, "/challenge-signature-correlations.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(abbas.cor.mat, file = paste0(figs.dir, "/abbas-signature-correlations.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cibersort.cor.mat, file = paste0(figs.dir, "/cibersort-signature-correlations.tsv"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


