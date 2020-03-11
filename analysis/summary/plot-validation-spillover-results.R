
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))

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
## validation.results.file <- "validation-results.csv"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path

res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
## print(head(subset(res, method == "cibersort" & subchallenge == "fine" & cell.type == "fibroblasts" & dataset.name == "DS5")))

## Limit to the final (of two) submissions for each group
## and the baseline methods (which were all submitted by Andrew L
## and hence only one of which is_latest)
## res <- subset(res, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))
## res <- subset(res, !is.na(measured))
## Unlike for the leaderboard data, validation data is only composed of measured cell types.
## So, if something is NA, that means it was not present.
flag <- is.na(res$measured)
res[flag, "measured"] <- 0

## Only keep the purified samples (these are in dataset 5)
flag <- res$dataset.name == "DS5"
res <- res[flag, ]

flag <- res$sample.id == "Breast"
res[flag, "sample.id"] <- "BRCA"

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
         "quantiseq" = "quanTIseq",
         "mcpcounter" = "MCP-counter",
         "epic" = "EPIC",
         "timer" = "TIMER",
         "xcell" = "xCell")
for(nm in names(trans)) {
    flag <- res$method == nm
    res[flag, "method"] <- trans[[nm]]
}

deconv.methods <- c("CIBERSORT", "EPIC", "quanTIseq")

flag <- grepl(res$sample.id, pattern="BM") | grepl(res$sample.id, pattern="RM")
if(any(flag)) {
  stop("Was not expecting andy BM or RM admixtures in dataset D5\n")
}

## Average replicates -- no so by removing _1 and _2
res$sample.id <- gsub(res$sample.id, pattern = "_1", replacement = "")
res$sample.id <- gsub(res$sample.id, pattern = "_2", replacement = "")

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
                   text = element_text(size = 20))
##                   title = element_text(size = 8))
    g <- g + ylab("Actual Cell Type") + xlab("Predicted Cell Type")
    if(normalized.score) {
        g <- g + scale_fill_gradient2("Normalized\nPrediction", 
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    } else {
        g <- g + scale_fill_gradient2("Prediction", limits = c(0,1),
                                      low = "red", high = "blue", mid = "white", na.value = "black")
    }
    g <- g + facet_wrap(~ method)
    g
}

png("validation-deconv-coarse-grained-spillover.png", width = 2 * 480)
g1 <- plot.cell.type.score.heatmap(subset(deconv.res, subchallenge == "coarse"))
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge (Validation)")
print(g1)
d <- dev.off()

png("validation-deconv-fine-grained-spillover.png", width = 2 * 480)
g2 <- plot.cell.type.score.heatmap(subset(deconv.res, subchallenge == "fine"))
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge (Validation)")
print(g2)
d <- dev.off()

non.deconv.res <- merge(non.deconv.res, purified.scores)
non.deconv.res$norm.score <- non.deconv.res$prediction / non.deconv.res$pure.score

png("validation-non-deconv-coarse-grained-spillover.png", width = 2 * 480)
g3 <- plot.cell.type.score.heatmap(subset(non.deconv.res, subchallenge == "coarse"),
                                  score.col = "norm.score", normalized.score = TRUE)
g3 <- g3 + ggtitle("Coarse-Grained Sub-Challenge (Validation)")
print(g3)
d <- dev.off()

png("validation-non-deconv-fine-grained-spillover.png", width = 2 * 480)
g4 <- plot.cell.type.score.heatmap(subset(non.deconv.res, subchallenge == "fine"),
                                  score.col = "norm.score", normalized.score = TRUE)
g4 <- g4 + ggtitle("Fine-Grained Sub-Challenge (Validation)")
print(g4)
d <- dev.off()

png("validation-coarse-grained-spillover.png", width = 2 * 480, height = 2 * 480)
g1 <- g1 + ggtitle("")
g3 <- g3 + ggtitle("")
title <- "Coarse-Grained Sub-Challenge (Validation)"
g <- grid.arrange(g1, g3, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
grid.draw(g)
d <- dev.off()

png("validation-fine-grained-spillover.png", width = 2 * 480, height = 2 * 480)
g2 <- g2 + ggtitle("")
g4 <- g4 + ggtitle("")
title <- "Fine-Grained Sub-Challenge (Validation)"
g <- grid.arrange(g2, g4, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
grid.draw(g)
d <- dev.off()

