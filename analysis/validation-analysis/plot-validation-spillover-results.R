
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(openxlsx))

## spillover -- cell type x to cell type y != x
##  - deconv (cibersort, epic, quantiseq); heatmap of fractions
##  - all (divide y by score of y for y-pure samples) --
##    i.e., score(actual.cell.type, predicted.cell.type), i.e., the score returned for predicted.cell.type
##          when only the actual.cell.type is present, should be normalized by
##          score(predicted.cell.type, predicted.cell.type)
##    NB: sometimes predicted.cell.type is a class that encompasses multiple actual.cell.types, in which
##        case we normalized by
##        mean_{cell.type in predicted.cell.type} score(cell.type, predicted.cell.type)
##        e.g., mean_{x in Memory_CD4_T_cells_1, Memory_CD4_T_cells_2, Naive_CD4_T_cells_1, Naive_CD4_T_cells_2, Tregs} score(x, CD4.T.cell)

## Get the validation results ("all_predictions.csv")
synLogin()
synId <- "syn21715094"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path
## validation.results.file <- "validation-results.csv"
## res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

prefix <- "all-genes-"
title.postfix <- "; All Genes"
## validation.results.file <- "all-predictions-all-genes.tsv"
csx.results.file <- "/Users/Brian/work/sage/deconvolution/Tumor-Deconvolution-Challenge/cibersortx/csx-all-gene-predictions.tsv"

synId <- "syn22331789"
csx.results.file <- synGet(synId, downloadFile = TRUE)$path

use.protein.coding.only <- FALSE
if(use.protein.coding.only) {

    prefix <- "only-pc-"
    title.postfix <- "; Protein-Coding Genes"
    title.postfix <- ""
    validation.results.file <- "all-predictions-pc-only.tsv"
    csx.results.file <- "/Users/Brian/work/sage/deconvolution/Tumor-Deconvolution-Challenge/cibersortx/csx-pc-only-predictions.tsv"
    
}
res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
csx.res <- read.table(csx.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
flag <- colnames(csx.res) == "method.name"
colnames(csx.res)[flag] <- "method"

measured <- unique(res[, c("subchallenge", "dataset.name", "sample.id", "cell.type", "measured")])
csx.res <- merge(csx.res, measured)

cols <- intersect(colnames(res), colnames(csx.res))

res <- rbind(res[, cols], csx.res[, cols])

synId <- "syn21590364"
path <- synGet(synId, downloadFile = TRUE)$path
val.coarse <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
val.coarse$subchallenge <- "coarse"

synId <- "syn21590365"
path <- synGet(synId, downloadFile = TRUE)$path
val.fine <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
val.fine$subchallenge <- "fine"

val.all <- rbind(val.coarse, val.fine)
flag <- val.all$sample.id == "Breast"
val.all[flag, "sample.id"] <- "BRCA"

## Only keep the purified samples (these are in dataset 5)
flag <- res$dataset.name == "DS5"
res <- res[flag, ]

flag <- res$sample.id == "Breast"
res[flag, "sample.id"] <- "BRCA"

res <- merge(res[, !(colnames(res) == "measured")], val.all)

## print(head(subset(res, method == "cibersort" & subchallenge == "fine" & cell.type == "fibroblasts" & dataset.name == "DS5")))

## Let's exclude memory.B.cells (which always have measured == 0, which causes problems with correlation)
res <- subset(res, !(cell.type == "memory.B.cells"))

## Limit to the final (of two) submissions for each group
## and the baseline methods (which were all submitted by Andrew L
## and hence only one of which is_latest)
## res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))
## res <- subset(res, !is.na(measured))
## Unlike for the leaderboard data, validation data is only composed of measured cell types.
## So, if something is NA, that means it was not present.
flag <- is.na(res$measured)
res[flag, "measured"] <- 0

average.replicates <- FALSE

sample.levels <- c(
    "Naive_B_cells",
    "Memory_CD4_T_cells",
    "Naive_CD4_T_cells",
    "Tregs",
    "Memory_CD8_T_cells",
    "Naive_CD8_T_cells",
    "NK_cells",
    "Neutrophils",
    "Monocytes",
    "Dendritic_cells",   
    "Macrophages",
    "Endothelial_cells",
    "Fibroblasts",
    "BRCA",
    "CRC"
)

if(average.replicates == FALSE) {
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
}

cell.type.levels <- c(
    "B.cells",
    "memory.B.cells",
    "naive.B.cells",
    "CD4.T.cells",
    "memory.CD4.T.cells",
    "naive.CD4.T.cells",
    "regulatory.T.cells",
    "CD8.T.cells",
    "memory.CD8.T.cells",
    "naive.CD8.T.cells",
    "NK.cells",
    "neutrophils",
    "monocytic.lineage",
    "monocytes",
    "myeloid.dendritic.cells",
    "macrophages",
    "endothelial.cells",
    "fibroblasts"
)

trans <-
    list("cibersort" = "CIBERSORT",
         "CIBERSORTx" = "CIBERSORTx",
         "quantiseq" = "quanTIseq",
         "mcpcounter" = "MCP-counter",
         "epic" = "EPIC",
         "timer" = "TIMER",
         "xcell" = "xCell")
for(nm in names(trans)) {
    flag <- res$method == nm
    res[flag, "method"] <- trans[[nm]]
}

deconv.methods <- c("CIBERSORT", "CIBERSORTx", "EPIC", "quanTIseq")

flag <- grepl(res$sample.id, pattern="BM") | grepl(res$sample.id, pattern="RM")
if(any(flag)) {
  stop("Was not expecting andy BM or RM admixtures in dataset D5\n")
}

## Average replicates -- do so by removing _1 and _2
if(average.replicates) {
    res$sample.id <- gsub(res$sample.id, pattern = "_1", replacement = "")
    res$sample.id <- gsub(res$sample.id, pattern = "_2", replacement = "")
}
res$sample.id <- factor(res$sample.id, levels = sample.levels)
    
res$cell.type <- factor(res$cell.type, levels = cell.type.levels)

res <-
    ddply(res,
          .variables = c("method", "subchallenge", "dataset.name", "sample.id", "cell.type"),
          .fun = function(df) {
              data.frame(prediction = mean(df$prediction),
                         measured = mean(df$measured))
          })

deconv.res <- subset(res, method %in% deconv.methods)
non.deconv.res <- subset(res, !(method %in% deconv.methods))

## Calculate the score for cell type y in y-purified cells
##    i.e., score(actual.cell.type, predicted.cell.type), i.e., the score returned for predicted.cell.type
##          when only the actual.cell.type is present, should be normalized by
##          score(predicted.cell.type, predicted.cell.type)
##    NB: sometimes predicted.cell.type is a class that encompasses multiple actual.cell.types, in which
##        case we normalized by
##        mean_{cell.type in predicted.cell.type} score(cell.type, predicted.cell.type)
##        e.g., mean_{x in Memory_CD4_T_cells_1, Memory_CD4_T_cells_2, Naive_CD4_T_cells_1, Naive_CD4_T_cells_2, Tregs} score(x, CD4.T.cell)
purified.scores <- subset(res, measured == 1)
purified.scores <-
    ddply(purified.scores,
          .variables = c("method", "subchallenge", "dataset.name", "cell.type"),
          .fun = function(df) { data.frame(pure.score = mean(df$prediction)) })

plot.cell.type.score.heatmap <- function(df, score.col = "prediction",
                                         normalized.score = FALSE) {
    g <- ggplot(data = df, aes_string(y = "sample.id", x = "cell.type", fill = score.col))
    g <- g + geom_tile()
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   text = element_text(size = 16))
##                   title = element_text(size = 8))
    g <- g + ylab("Purified Sample") + xlab("Predicted Cell Type")
    if(normalized.score) {
        g <- g + scale_fill_gradient2("Normalized\nPrediction", 
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    } else {
        g <- g + scale_fill_gradient2("Prediction", limits = c(0,1),
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    }
    g <- g + facet_wrap(~ method, nrow = 1)
    g
}

all.res <- merge(res, purified.scores, all.x = TRUE)
all.res$norm.score <- all.res$prediction / all.res$pure.score

coarse.res <- subset(all.res, subchallenge == "coarse")
fine.res <- subset(all.res, subchallenge == "fine")

non.deconv.methods <- unique(as.character(all.res$method))
non.deconv.methods <- non.deconv.methods[!(non.deconv.methods %in% deconv.methods)]
names(deconv.methods) <- deconv.methods
names(non.deconv.methods) <- non.deconv.methods
l_ply(deconv.methods,
      .fun = function(meth) {
          g <- plot.cell.type.score.heatmap(subset(coarse.res, method == meth))
          g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\n(Validation", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
          file <- paste0(prefix, "validation-spillover-coarse-grained-", meth, ".png")
          png(file)
          print(g)
          d <- dev.off()
          
          g <- plot.cell.type.score.heatmap(subset(fine.res, method == meth))
          g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\n(Validation", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
          file <- paste0(prefix, "validation-spillover-fine-grained-", meth, ".png")
          png(file)
          print(g)
          d <- dev.off()
      })

l_ply(non.deconv.methods,
      .fun = function(meth) {
          g <- plot.cell.type.score.heatmap(subset(coarse.res, method == meth),
                                            score.col = "norm.score", normalized.score = TRUE)
          g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\n(Validation", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
          file <- paste0(prefix, "validation-spillover-coarse-grained-", meth, ".png")
          png(file)
          print(g)
          d <- dev.off()
          
          g <- plot.cell.type.score.heatmap(subset(fine.res, method == meth),
                                            score.col = "norm.score", normalized.score = TRUE)
          g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\n(Validation", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
          file <- paste0(prefix, "validation-spillover-fine-grained-", meth, ".png")
          png(file)
          print(g)
          d <- dev.off()
      })


png(paste0(prefix, "validation-spillover-deconv-coarse-grained.png"), width = 2 * 480)
g1 <- plot.cell.type.score.heatmap(subset(deconv.res, subchallenge == "coarse"))
g1 <- g1 + ggtitle(paste0("Coarse-Grained Sub-Challenge (Validation", title.postfix, ")"))
print(g1)
d <- dev.off()

png(paste0(prefix, "validation-spillover-deconv-fine-grained.png"), width = 2 * 480)
g2 <- plot.cell.type.score.heatmap(subset(deconv.res, subchallenge == "fine"))
g2 <- g2 + ggtitle(paste0("Fine-Grained Sub-Challenge (Validation", title.postfix, ")"))
print(g2)
d <- dev.off()

non.deconv.res <- merge(non.deconv.res, purified.scores)
non.deconv.res$norm.score <- non.deconv.res$prediction / non.deconv.res$pure.score

png(paste0(prefix, "validation-spillover-non-deconv-coarse-grained.png"), width = 2 * 480)
g3 <- plot.cell.type.score.heatmap(subset(non.deconv.res, subchallenge == "coarse"),
                                  score.col = "norm.score", normalized.score = TRUE)
g3 <- g3 + ggtitle(paste0("Coarse-Grained Sub-Challenge (Validation", title.postfix, ")"))
print(g3)
d <- dev.off()

png(paste0(prefix, "validation-spillover-non-deconv-fine-grained.png"), width = 2 * 480)
g4 <- plot.cell.type.score.heatmap(subset(non.deconv.res, subchallenge == "fine"),
                                  score.col = "norm.score", normalized.score = TRUE)
g4 <- g4 + ggtitle(paste0("Fine-Grained Sub-Challenge (Validation", title.postfix, ")"))
print(g4)
d <- dev.off()

png(paste0(prefix, "validation-spillover-all-coarse-grained.png"), width = 4 * 480, height = 2 * 480)
g1 <- g1 + ggtitle("")
g3 <- g3 + ggtitle("")
g1 <- g1 + theme(axis.text.x = element_text(size = 15))
g3 <- g3 + theme(axis.text.x = element_text(size = 15))
title <- paste0("Coarse-Grained Sub-Challenge (Validation", title.postfix, ")")
g <- grid.arrange(g1, g3, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
grid.draw(g)
d <- dev.off()

png(paste0(prefix, "validation-spillover-all-fine-grained.png"), width = 4 * 480, height = 2 * 480)
g2 <- g2 + ggtitle("")
g4 <- g4 + ggtitle("")
g2 <- g2 + theme(axis.text.x = element_text(size = 15))
g4 <- g4 + theme(axis.text.x = element_text(size = 15))
title <- paste0("Fine-Grained Sub-Challenge (Validation", title.postfix, ")")
g <- grid.arrange(g2, g4, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
grid.draw(g)
d <- dev.off()

