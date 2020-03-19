
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ComplexHeatmap))
suppressPackageStartupMessages(p_load(openxlsx))

synLogin()

## Get the _transcript_ expression matrix
synId <- "syn21574261"
trans.expr.mat.file <- synGet(synId, downloadFile = TRUE)$path
trans.expr.mat <- read.table(trans.expr.mat.file, sep="\t", header=TRUE)
cd45ro.enst.id <- "ENST00000348564.11"
rownames(trans.expr.mat) <- as.character(trans.expr.mat$transcript)
trans.expr.mat <- trans.expr.mat[, !(colnames(trans.expr.mat) == "transcript")]

## Get the expression matrix
synId <- "syn21576632"
expr.mat.file <- synGet(synId, downloadFile = TRUE)$path

## CD45RO = ENST00000348564.11 = PTPRC-201

expr.mat <- read.table(expr.mat.file, sep=",", header=TRUE)
rownames(expr.mat) <- as.character(expr.mat$Gene)
expr.mat <- expr.mat[, !(colnames(expr.mat) == "Gene")]

## Only keep the purified samples
colnames(expr.mat)[colnames(expr.mat) == "Breast"] <- "BRCA"
colnames(trans.expr.mat)[colnames(trans.expr.mat) == "Breast"] <- "BRCA"
cols <- colnames(expr.mat)
##cols <- cols[!(grepl(cols, pattern="BM"))]
##cols <- cols[!(grepl(cols, pattern="RM"))]


sample.levels <- c(
    "Naive_B_cells_1",
    "Memory_CD4_T_cells_1",
    "Memory_CD4_T_cells_2",        
    "Naive_CD4_T_cells_1",
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
trans.expr.mat <- trans.expr.mat[, cols]

expr.mat <- rbind(expr.mat, CD45RO = as.numeric(trans.expr.mat[cd45ro.enst.id,]))

## Read in the marker table
marker.file <- "immune-cell-markers.xlsx"
marker.tbl <- read.xlsx(marker.file, sheet=1)

plot.marker.heatmap <- function(mat, marker.tbl, marker.gene.id.col = "marker.symbol", zscore = TRUE) {
    if(zscore) {
        mat <- t(scale(t(mat), scale = TRUE, center = TRUE))
    }
    rownames(marker.tbl) <- as.character(marker.tbl[, marker.gene.id.col])
    common <- intersect(rownames(mat), as.character(marker.tbl[, marker.gene.id.col]))
    mtbl <- marker.tbl[rownames(marker.tbl) %in% common,]
    mat <- mat[as.character(mtbl[, marker.gene.id.col]),]
    
    split <- factor(mtbl$cell.type, levels = unique(mtbl$cell.type))
    h <- Heatmap(mat, cluster_rows = FALSE, show_column_names = TRUE, split = split, cluster_columns = FALSE,
                 row_title_rot = 0, cluster_row_slices = FALSE, column_names_rot = 90)
    h
}

h1 <- plot.marker.heatmap(as.matrix(expr.mat), marker.tbl)
h1@column_title <- "All genes"
##draw(h1, column_title = "All genes")

source("../utils.R")
pc.expr.mat <- limit.matrix.to.protein.coding(expr.mat)
pc.expr.mat <- rbind(pc.expr.mat, CD45RO = as.numeric(trans.expr.mat[cd45ro.enst.id,]))

h2 <- plot.marker.heatmap(as.matrix(pc.expr.mat), marker.tbl)
##draw(h1, column_title = "All genes")
h2@column_title <- "Protein-coding genes"

png("purified-samples-marker-heatmap.png", width = 2 * 480)
draw(h1 + h2)
d <- dev.off()

png("purified-samples-marker-heatmap-all-genes.png")
draw(h1)
d <- dev.off()

png("purified-samples-marker-heatmap-protein-coding-genes.png")
draw(h2)
d <- dev.off()

          
