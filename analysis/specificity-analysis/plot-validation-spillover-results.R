
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

## Get the sensitivity results
synLogin()
synId <- "syn22491951"
validation.results.file <- synGet(synId, downloadFile = TRUE)$path

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

res <- read.table(validation.results.file, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)

source("../utils.R")

res <- assign.result.team.names.and.rounds(res, error.fun = warning)

include.csx.results <- FALSE
if(include.csx.results) {
  synId <- "syn22331789"
  csx.results.file <- synGet(synId, downloadFile = TRUE)$path
  csx.res <- read.table(csx.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
  flag <- colnames(csx.res) == "method.name"
  colnames(csx.res)[flag] <- "method"

  measured <- unique(res[, c("subchallenge", "dataset.name", "sample.id", "cell.type", "measured")])
  csx.res <- merge(csx.res, measured)
  
  cols <- intersect(colnames(res), colnames(csx.res))

  res <- rbind(res[, cols], csx.res[, cols])
}

flag <- res$sample.id == "Breast"
res[flag, "sample.id"] <- "BRCA"

## Let's exclude memory.B.cells (which always have measured == 0, which causes problems with correlation)
res <- subset(res, !(cell.type == "memory.B.cells"))

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

deconv.methods <- c("CIBERSORT", "EPIC", "quanTIseq")
if(include.csx.results) {
    deconv.methods <- c("CIBERSORTx", deconv.methods)
}

flag <- grepl(res$sample.id, pattern="BM") | grepl(res$sample.id, pattern="RM")
if(any(flag)) {
  stop("Was not expecting any BM or RM admixtures in dataset D5\n")
}

## Average replicates -- do so by removing _1 and _2
if(average.replicates) {
    res$sample.id <- gsub(res$sample.id, pattern = "_1", replacement = "")
    res$sample.id <- gsub(res$sample.id, pattern = "_2", replacement = "")
}
res$sample.id <- factor(res$sample.id, levels = sample.levels)
    
res$cell.type <- factor(res$cell.type, levels = cell.type.levels)

if(FALSE) {
    ## This appears to be averaging over submission round
    res <-
        ddply(res,
              .variables = c("method.name", "subchallenge", "dataset.name", "sample.id", "cell.type"),
              .fun = function(df) {
                  data.frame(prediction = mean(df$prediction),
                             measured = mean(df$measured))
              })
} # if(FALSE)

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
    g <- g + facet_wrap(~ method.name, nrow = 1)
    g
}

plot.spillover.results <- function(res.input,
                                   method.name.col, 
                                   subchallenge.col, 
                                   round.col, round = "latest",
                                   postfix) {

    submitter.tbl <- unique(res.input[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
    or <- order(submitter.tbl[, round.col])
    submitter.tbl <- submitter.tbl[or, ]
    submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
    flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
    submitter.tbl <- submitter.tbl[flag, ]
    flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
    submitter.tbl <- submitter.tbl[flag, ]

    res <- merge(res.input, submitter.tbl, by = c(method.name.col, subchallenge.col, round.col))

    round.text <- ""
    if(round == "latest") {
        round.text <- "Latest Round"
    } else if (round == "1") {
        round.text <- "Round 1"
    } else {
        round.text <- paste0("Latest Round up to Round ", round)
    }

    title.postfix <- round.text
    
    flag <- res[, method.name.col] %in% deconv.methods
    deconv.res <- res[flag, ]
    non.deconv.res <- res[!flag, ]
    
    ## Calculate the score for cell type y in y-purified cells
    ## OLD APPROACH:
    ##    i.e., score(actual.cell.type, predicted.cell.type), i.e., the score returned for predicted.cell.type
    ##          when only the actual.cell.type is present, should be normalized by
    ##          score(predicted.cell.type, predicted.cell.type)
    ##    NB: sometimes predicted.cell.type is a class that encompasses multiple actual.cell.types, in which
    ##        case we normalized by
    ##        mean_{cell.type in predicted.cell.type} score(cell.type, predicted.cell.type)
    ##        e.g., mean_{x in Memory_CD4_T_cells_1, Memory_CD4_T_cells_2, Naive_CD4_T_cells_1, Naive_CD4_T_cells_2, Tregs} score(x, CD4.T.cell)
    ## NEW APPROACH:
    ## normalized.score(actual, predicted) = [ score(actual, predicted) - min_predicted' score(actual = 0, predicted') ] /
    ##                                       [ max_predicted' score(actual = 1, predicted') - min_predicted' score(actual = 0, predicted') ]
    ## where
    ## min_predicted' score(actual = 0, predicted') is effectively the background for the cell type with actual = 1
    ## NEWEST APPROACH:
    ## don't limit max_predicted' above to cases where actual = 1.
    ## e.g., TIMER returns a higher score for macrophages when there are no macrophages in the sample, as opposed to
    ## when it is pure macrophages
    normalizing.scores <- subset(res, measured == 1)
    normalizing.scores <-
        ddply(normalizing.scores,
              .variables = c(method.name.col, subchallenge.col, "dataset.name", "cell.type"),
              .fun = function(df) { data.frame(pure.score = mean(df$prediction), 
                                               max.purified.score = max(df$prediction)) })

    tmp <-
        ddply(res,
              .variables = c(method.name.col, subchallenge.col, "dataset.name", "cell.type"),
              .fun = function(df) { data.frame(max.score = max(df$prediction)) })
    

    normalizing.scores <- merge(normalizing.scores, tmp)

    background.scores <- subset(res, measured == 0)
    background.scores <-
        ddply(background.scores,
              .variables = c(method.name.col, subchallenge.col, "dataset.name", "cell.type"),
              .fun = function(df) { data.frame(min.background.score = min(df$prediction)) })
    
    normalizing.scores <- merge(normalizing.scores, background.scores)
    
    all.res <- merge(res, normalizing.scores, all.x = TRUE)
    ## all.res$norm.score <- all.res$prediction / all.res$pure.score
    ## all.res$norm.score <- ( all.res$prediction - all.res$min.background.score ) / ( all.res$max.purified.score - all.res$min.background.score )
    all.res$norm.score <- ( all.res$prediction - all.res$min.background.score ) / ( all.res$max.score - all.res$min.background.score )

    flag <- all.res[, subchallenge.col] == "coarse"
    coarse.res <- all.res[flag, ]
    flag <- all.res[, subchallenge.col] == "fine"
    fine.res <- all.res[flag, ]
    
    non.deconv.methods <- unique(as.character(all.res[, method.name.col]))
    non.deconv.methods <- non.deconv.methods[!(non.deconv.methods %in% deconv.methods)]
    
    ## non.deconv.methods <- c("xCell", "MCP-counter", "TIMER")
    
    names(deconv.methods) <- deconv.methods
    names(non.deconv.methods) <- non.deconv.methods
    l_ply(deconv.methods,
          .fun = function(meth) {
              flag <- coarse.res[, method.name.col] == meth
              g <- plot.cell.type.score.heatmap(coarse.res[flag, ])
              g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\n(", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
              file <- paste0("spillover-coarse-grained-", make.names(meth), postfix, ".png")
              png(file)
              print(g)
              d <- dev.off()

              flag <- fine.res[, method.name.col] == meth              
              g <- plot.cell.type.score.heatmap(fine.res[flag,])
              g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\n(", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
              file <- paste0("spillover-fine-grained-", make.names(meth), postfix, ".png")
              png(file)
              print(g)
              d <- dev.off()
          })
    
    l_ply(non.deconv.methods,
          .fun = function(meth) {
              flag <- coarse.res[, method.name.col] == meth                            
              tmp <- coarse.res[flag, ]
              if(nrow(tmp) > 0) {
                  g <- plot.cell.type.score.heatmap(tmp, score.col = "norm.score", normalized.score = TRUE)
                  g <- g + ggtitle(paste0("Coarse-Grained Sub-Challenge\n(", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
                  file <- paste0("spillover-coarse-grained-", make.names(meth), postfix, ".png")
                  png(file)
                  print(g)
                  d <- dev.off()
              }

              flag <- fine.res[, method.name.col] == meth                                          
              tmp <- fine.res[flag, ]
              if(nrow(tmp) > 0) {          
                  g <- plot.cell.type.score.heatmap(tmp, score.col = "norm.score", normalized.score = TRUE)
                  g <- g + ggtitle(paste0("Fine-Grained Sub-Challenge\n(", title.postfix, ")")) + theme(plot.title = element_text(hjust = 0.5))
                  file <- paste0("spillover-fine-grained-", make.names(meth), postfix, ".png")
                  png(file)
                  print(g)
                  d <- dev.off()
              }
          })
    
    
    png(paste0("spillover-deconv-coarse-grained", postfix, ".png"), width = 2 * 480)
    flag <- deconv.res[, subchallenge.col] == "coarse"
    g1 <- plot.cell.type.score.heatmap(deconv.res[flag, ])
    g1 <- g1 + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g1)
    d <- dev.off()
    
    png(paste0("spillover-deconv-fine-grained", postfix, ".png"), width = 2 * 480)
    flag <- deconv.res[, subchallenge.col] == "fine"    
    g2 <- plot.cell.type.score.heatmap(deconv.res[flag, ])
    g2 <- g2 + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g2)
    d <- dev.off()

    non.deconv.res <- merge(non.deconv.res, normalizing.scores)
    ## non.deconv.res$norm.score <- non.deconv.res$prediction / non.deconv.res$pure.score
    ## non.deconv.res$norm.score <- ( non.deconv.res$prediction - non.deconv.res$min.background.score ) / ( non.deconv.res$max.purified.score - non.deconv.res$min.background.score )
    non.deconv.res$norm.score <- ( non.deconv.res$prediction - non.deconv.res$min.background.score ) / ( non.deconv.res$max.score - non.deconv.res$min.background.score )
    flag <- (non.deconv.res$prediction == 0) & (non.deconv.res$max.purified.score == 0) & (non.deconv.res$min.background.score == 0)
    non.deconv.res[flag, "norm.score"] <- 0
    
    
    png(paste0("spillover-non-deconv-coarse-grained", postfix, ".png"), width = 2 * 480)
    flag <- non.deconv.res[, subchallenge.col] == "coarse"        
    sub <- non.deconv.res[flag, ]
    g3 <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE)
    g3 <- g3 + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g3)
    d <- dev.off()
    
    png(paste0("spillover-non-deconv-fine-grained", postfix, ".png"), width = 2 * 480)
    flag <- non.deconv.res[, subchallenge.col] == "fine"        
    sub <- non.deconv.res[flag, ]
    g4 <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE)
    g4 <- g4 + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g4)
    d <- dev.off()
    
    png(paste0("spillover-all-coarse-grained", postfix, ".png"), width = 4 * 480, height = 2 * 480)
    g1 <- g1 + ggtitle("")
    g3 <- g3 + ggtitle("")
    g1 <- g1 + theme(axis.text.x = element_text(size = 15))
    g3 <- g3 + theme(axis.text.x = element_text(size = 15))
    title <- paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")")
    g <- grid.arrange(g1, g3, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
    grid.draw(g)
    d <- dev.off()
    
    png(paste0("spillover-all-fine-grained", postfix, ".png"), width = 4 * 480, height = 2 * 480)
    g2 <- g2 + ggtitle("")
    g4 <- g4 + ggtitle("")
    g2 <- g2 + theme(axis.text.x = element_text(size = 15))
    g4 <- g4 + theme(axis.text.x = element_text(size = 15))
    title <- paste0("Fine-Grained Sub-Challenge (", title.postfix, ")")
    g <- grid.arrange(g2, g4, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
    grid.draw(g)
    d <- dev.off()

}

results <- list()
for(round in c("1", "2", "3", "latest")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    plot.spillover.results(res, method.name.col, subchallenge.col, round.col = "submission", round = round, postfix = postfix)
}
