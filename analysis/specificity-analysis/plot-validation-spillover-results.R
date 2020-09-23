
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(openxlsx))
suppressPackageStartupMessages(p_load(reshape2))

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

include.csx.results <- TRUE
if(include.csx.results) {
    synIds <- list("coarse" = "syn22725972", "fine" = "syn22726111")
    nms <- names(synIds)
    names(nms) <- nms
    csx.res <-
        ldply(nms,
              .fun = function(nm) {
                  synId <- synIds[[nm]]
                  csx.results.file <- synGet(synId, downloadFile = TRUE)$path
                  tbl <- read.table(csx.results.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors = FALSE)
                  flag <- tbl$subchallenge == nm
                  tbl <- tbl[flag, ]
                  tbl$objectId <- "CIBERSORTx"
                  tbl$repo_name <- "CIBERSORTx"
                  tbl$comparator <- TRUE
                  tbl$submitterId <- NA
                  tbl1 <- tbl
                  tbl2 <- tbl
                  tbl1$submission <- "1"
                  tbl2$submission <- "latest"                  
                  rbind(tbl1, tbl2, stringsAsFactors = FALSE)
              })
    measured <- unique(res[, c("subchallenge", "dataset.name", "sample.id", "cell.type", "measured")])
    csx.res <- merge(csx.res, measured)
  
    cols <- intersect(colnames(res), colnames(csx.res))

    res <- rbind(res[, cols], csx.res[, cols])
}


flag <- res$sample.id == "Breast"
res[flag, "sample.id"] <- "BRCA"

## Let's exclude memory.B.cells (which always have measured == 0, which causes problems with correlation)
res <- subset(res, !(cell.type == "memory.B.cells"))

res[, cell.type.col] <- as.character(res[, cell.type.col])
res <- rename.cell.types(res, from.col = cell.type.col, to.col = cell.type.col)


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

df <- data.frame(cell.type = cell.type.levels, stringsAsFactors = FALSE)
df <- rename.cell.types(df, from.col = "cell.type", to.col = "cell.type")
cell.type.levels <- as.character(df$cell.type)


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

## Order the samples and the cell.type.levels consistently (so that we get a diagonal
## in the spillover heatmaps)
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
                                         normalized.score = FALSE, nrow = 1, method.name.col = "method.name") {

    ## This re-ordering o works with as.table = FALSE to respect the
    ## the ordering we want from cell.type.levels
    if(nrow > 1) {
        method.levels <- sort(unique(df[, method.name.col]), decreasing = TRUE)
        ntot <- length(method.levels)
        ncol <- ceiling(ntot / nrow)
        o <- unlist(llply(1:(nrow-1), .fun = function(i) (i*ncol):(1+((i-1)*ncol))))
        o <- c(o, ntot:(max(o)+1))
        df[, method.name.col] <- factor(df[, method.name.col], levels = method.levels[o])
    }
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
    g <- g + facet_wrap(method.name.col, as.table = FALSE, nrow = nrow)
    g
}

perform.spillover.analysis <- function(res.input,
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

    flag <- (all.res$prediction == 0) & (all.res$max.purified.score == 0) & (all.res$min.background.score == 0)
    all.res[flag, "norm.score"] <- 0

    ## The teams we invited to participate in the collaborative phase
    top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
    priority.methods <- unique(c(top.performers, unique(subset(all.res, comparator==TRUE)[, method.name.col])))
    
    flag <- all.res[, method.name.col] %in% deconv.methods
    deconv.res <- all.res[flag, ]
    non.deconv.res <- all.res[!flag, ]
    
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
    

    l <- list("coarse" = subset(coarse.res, measured == 0), "fine" = subset(fine.res, measured == 0))
    mean.norm.score.by.cell.type.method <-
            llply(l,
                  .fun = function(df) {
                      na.rm <- FALSE
                      
                      ret <- ddply(df, .variables = c(method.name.col, cell.type.col),
                                   .fun = function(df) {
                                       data.frame(norm.score = mean(df[, "norm.score"], na.rm=na.rm))
                                   })

                      ## Introduce NAs for missing entries
                      ret <- melt(acast(ret, as.formula(paste0(cell.type.col, "~", method.name.col))))
                      colnames(ret) <- c(cell.type.col, method.name.col, "norm.score")
                      ret
                  })

    
    png(paste0("spillover-deconv-coarse-grained", postfix, ".png"), width = 2 * 480)
    flag <- deconv.res[, subchallenge.col] == "coarse"
    g1 <- plot.cell.type.score.heatmap(deconv.res[flag, ])
    g1 <- g1 + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g1)
    d <- dev.off()

    g.strip.deconv.coarse <- plot.strip.plots(subset(deconv.res[flag, ], measured == 0),
                                          id.var = method.name.col, cell.type.var = cell.type.col,
                                          var = "prediction", label = "Prediction", col.summary.fun = "mean",
                                          order.decreasing = TRUE)
    g.strip.deconv.coarse <- g.strip.deconv.coarse + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    
    png(paste0("spillover-deconv-fine-grained", postfix, ".png"), width = 2 * 480)
    flag <- deconv.res[, subchallenge.col] == "fine"    
    g2 <- plot.cell.type.score.heatmap(deconv.res[flag, ])
    g2 <- g2 + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g2)
    d <- dev.off()

    g.strip.deconv.fine <- plot.strip.plots(subset(deconv.res[flag, ], measured == 0),
                                            id.var = method.name.col, cell.type.var = cell.type.col,
                                            var = "prediction", label = "Prediction", col.summary.fun = "mean",
                                            order.decreasing = TRUE)
    g.strip.deconv.fine <- g.strip.deconv.fine + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    
    png(paste0("spillover-non-deconv-coarse-grained", postfix, ".png"), width = 2 * 480)
    flag <- non.deconv.res[, subchallenge.col] == "coarse"        
    sub <- non.deconv.res[flag, ]
    g3 <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE)
    g3 <- g3 + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g3)
    d <- dev.off()

    g.strip.non.deconv.coarse <- plot.strip.plots(subset(sub, measured == 0),
                                                  id.var = method.name.col, cell.type.var = cell.type.col,
                                                  var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                                  order.decreasing = TRUE)
    g.strip.non.deconv.coarse <- g.strip.non.deconv.coarse + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))

    
    png(paste0("spillover-non-deconv-fine-grained", postfix, ".png"), width = 2 * 480)
    flag <- non.deconv.res[, subchallenge.col] == "fine"        
    sub <- non.deconv.res[flag, ]
    g4 <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE)
    g4 <- g4 + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g4)
    d <- dev.off()

    g.strip.non.deconv.fine <- plot.strip.plots(subset(sub, measured == 0),
                                                id.var = method.name.col, cell.type.var = cell.type.col,
                                                var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                                order.decreasing = TRUE)
    g.strip.non.deconv.fine <- g.strip.deconv.fine + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))


    
    png(paste0("spillover-all-coarse-grained", postfix, ".png"), width = 4 * 480, height = 2 * 480)
    g1 <- g1 + ggtitle("")
    g3 <- g3 + ggtitle("")
    g1 <- g1 + theme(axis.text.x = element_text(size = 15))
    g3 <- g3 + theme(axis.text.x = element_text(size = 15))
    title <- paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")")
    g <- grid.arrange(g1, g3, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
    grid.draw(g)
    d <- dev.off()

    png(paste0("spillover-all-scores-coarse-grained", postfix, ".png"), width = 2 * 480)
    sub <- coarse.res
    g.all.coarse <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE, nrow = 2)
    g.all.coarse <- g.all.coarse + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.coarse)
    d <- dev.off()

    ## Top-performers only
    png(paste0("spillover-all-scores-coarse-grained-top", postfix, ".png"), width = 2 * 480)
    sub <- coarse.res
    flag <- sub[, method.name.col] %in% priority.methods
    sub <- sub[flag, ]
    g.all.coarse.top <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE, nrow = 2)
    g.all.coarse.top <- g.all.coarse.top + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.coarse.top)
    d <- dev.off()
    
    png(paste0("spillover-all-coarse-grained-strip", postfix, ".png"), width = 4 * 480, height = 2 * 480)
    g.strip.deconv.coarse <- g.strip.deconv.coarse + ggtitle("")
    g.strip.non.deconv.coarse <- g.strip.non.deconv.coarse + ggtitle("")
    g.strip.deconv.coarse <- g.strip.deconv.coarse + theme(axis.text.x = element_text(size = 15))
    g.strip.non.deconv.coarse <- g.strip.non.deconv.coarse + theme(axis.text.x = element_text(size = 15))
    title <- paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")")
    g <- grid.arrange(g.strip.deconv.coarse, g.strip.non.deconv.coarse, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
    grid.draw(g)
    d <- dev.off()

    png(paste0("spillover-all-scores-coarse-grained-strip", postfix, ".png"), width = 2 * 480)
    sub <- coarse.res
    g.all.strip.coarse <- plot.strip.plots(subset(sub, measured == 0),
                                           id.var = method.name.col, cell.type.var = cell.type.col,
                                           var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                           order.decreasing = TRUE)
    g.all.strip.coarse <- g.all.strip.coarse + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.strip.coarse)
    d <- dev.off()

    ## Top-performers only
    png(paste0("spillover-all-scores-coarse-grained-strip-top", postfix, ".png"), width = 2 * 480)
    sub <- coarse.res
    flag <- sub[, method.name.col] %in% priority.methods
    sub <- sub[flag, ]
    g.all.strip.coarse.top <- plot.strip.plots(subset(sub, measured == 0),
                                               id.var = method.name.col, cell.type.var = cell.type.col,
                                               var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                               order.decreasing = TRUE)
    g.all.strip.coarse.top <- g.all.strip.coarse.top + ggtitle(paste0("Coarse-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.strip.coarse.top)
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

    png(paste0("spillover-all-scores-fine-grained", postfix, ".png"), width = 2 * 480)
    sub <- fine.res
    g.all.fine <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE, nrow = 2)
    g.all.fine <- g.all.fine + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.fine)
    d <- dev.off()

    ## Top-performers only    
    png(paste0("spillover-all-scores-fine-grained-top", postfix, ".png"), width = 2 * 480)
    sub <- fine.res
    flag <- sub[, method.name.col] %in% priority.methods
    sub <- sub[flag, ]
    g.all.fine.top <- plot.cell.type.score.heatmap(sub, score.col = "norm.score", normalized.score = TRUE, nrow = 2)
    g.all.fine.top <- g.all.fine.top + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.fine.top)
    d <- dev.off()

    png(paste0("spillover-all-fine-grained-strip", postfix, ".png"), width = 4 * 480, height = 2 * 480)
    g.strip.deconv.fine <- g.strip.deconv.fine + ggtitle("")
    g.strip.non.deconv.fine <- g.strip.non.deconv.fine + ggtitle("")
    g.strip.deconv.fine <- g.strip.deconv.fine + theme(axis.text.x = element_text(size = 15))
    g.strip.non.deconv.fine <- g.strip.non.deconv.fine + theme(axis.text.x = element_text(size = 15))
    title <- paste0("Fine-Grained Sub-Challenge (", title.postfix, ")")
    g <- grid.arrange(g.strip.deconv.fine, g.strip.non.deconv.fine, nrow = 2, top = textGrob(title, gp = gpar(fontsize = 25)))
    grid.draw(g)
    d <- dev.off()

    png(paste0("spillover-all-scores-fine-grained-strip", postfix, ".png"), width = 2 * 480)
    sub <- fine.res
    g.all.strip.fine <- plot.strip.plots(subset(sub, measured == 0),
                                         id.var = method.name.col, cell.type.var = cell.type.col,
                                         var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                         order.decreasing = TRUE)
    g.all.strip.fine <- g.all.strip.fine + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.strip.fine)
    d <- dev.off()

    ## Top-performers only
    png(paste0("spillover-all-scores-fine-grained-strip-top", postfix, ".png"), width = 2 * 480)
    sub <- fine.res
    flag <- sub[, method.name.col] %in% priority.methods
    sub <- sub[flag, ]
    g.all.strip.fine.top <- plot.strip.plots(subset(sub, measured == 0),
                                             id.var = method.name.col, cell.type.var = cell.type.col,
                                             var = "norm.score", label = "Normalized Score", col.summary.fun = "mean",
                                             order.decreasing = TRUE)
    g.all.strip.fine.top <- g.all.strip.fine.top + ggtitle(paste0("Fine-Grained Sub-Challenge (", title.postfix, ")"))
    print(g.all.strip.fine.top)
    d <- dev.off()
    
    for(sub.challenge in names(mean.norm.score.by.cell.type.method)) {
        means <- mean.norm.score.by.cell.type.method[[sub.challenge]]
        g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "norm.score",
                                                col.summary.fun = "min", cor.type.label = "Normalized\nScore",
                                                limits = c(0, 1),
                                                order.decreasing = TRUE)
        png(paste0("spillover-all-scores-", sub.challenge, "-grained-heatmap", postfix, ".png"), width = 2 * 480)
        print(g)
        d <- dev.off()
    }
    
    return(all.res)
}

results <- list()
for(round in c("1", "2", "3", "latest")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    results[[round]] <- perform.spillover.analysis(res, method.name.col, subchallenge.col,
                                                   round.col = "submission", round = round, postfix = postfix)
    
}


