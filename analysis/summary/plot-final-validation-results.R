suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ComplexHeatmap))
suppressPackageStartupMessages(p_load(reshape2))

source("../utils.R")

synLogin()

## Read the leaderboards
sub.challenges <- list("coarse" = "coarse", "fine" = "fine")
synIds <- list("coarse" = "syn21822681", "fine" = "syn21822682")
lb.dfs <- llply(synIds,
                .fun = function(synId) {
                    synTableQuery(paste0("SELECT * FROM ", synId))$asDataFrame()
                })

submitter.id.col <- "submitterId"
object.id.col <- "objectId"
method.name.col <- "repo_name"

lb.dfs <- llply(lb.dfs,
                .fun = function(df) {
                    df$modelId <- paste0(as.character(df$repo_name), "-",
                                         as.character(df$submitter))
                    df <- simplify.submitter.names(df)
                    df <- assign.baseline.names(df, from.col = "repo_name", to.col = "submitter")
                    ## For some reason cibersortX has 0's for mem and naive CD8s, but it should be NA
                    ## (it doesn't call these)
                    flag <- (df$submitter == "CIBERSORTx") & (df$celltype %in% c("naive.CD8.T.cells", "memory.CD8.T.cells"))
                    df[flag, "metric_value"] <- NA
                    context.cols <- c(submitter.id.col)
                    df <- assign.submission.rounds(df, object.id.col, context.cols, method.name.col,
                                                   assign.latest = TRUE)
                    df
                })

for(nm in names(synIds)) {
    write.table(file = paste0("deconvolution-final-validation-", nm, ".tsv"), lb.dfs[[nm]], row.names = FALSE, col.names = TRUE,
                sep = "\t", quote = FALSE)
}

## NB: pearson is primary metric and spearman is secondary

## Plot the grand mean
plot.grand.mean <- function(df, col = "dataset") {
    df <- df[df[,col] == "Grand mean", ]
    df <- df[!is.na(df$metric_value),]
    
    ## Order the methods by pearson correlation (the primary metric)
    df.pearson <- subset(df, metric == "pearson")
    df.pearson <- df.pearson[order(df.pearson$metric_value, decreasing=TRUE),]
    ## Temporary workaround for fact that there are multiple is_latest = TRUE
    ## flag <- !duplicated(df.pearson$submitter)
    flag <- df.pearson$submission == "latest"
    df.pearson <- df.pearson[flag,]
    
    lvls <- unique(df.pearson$submitter)

    ## Just in case we've missed some methods (within an undefined pearson correlation!?)
    ## add them, as ordered by spearman correlation
    df.spearman <- subset(df, metric == "spearman")
    df.spearman <- df.spearman[order(df.spearman$metric_value, decreasing=TRUE),]
    ## Temporary workaround for fact that there are multiple is_latest = TRUE
    ## flag <- !duplicated(df.spearman$submitter)
    flag <- df.spearman$submission == "latest"
    df.spearman <- df.spearman[flag,]
    spearman.lvls <- unique(df.spearman$submitter)    
    lvls <- rev(c(lvls, spearman.lvls[!(spearman.lvls %in% lvls)]))

    df$submitter <- factor(df$submitter, levels = lvls)
    metrics <- list("pearson" = "pearson", "spearman" = "spearman")

    gs <-
        llply(metrics,
              .fun = function(met) {
                  df.met <- subset(df, metric == met)
                  ## latest submission will duplicate one of submission 1, 2, or 3
                  df.met <- subset(df.met, submission != "latest")
                  df.met <- df.met[order(df.met$metric_value, decreasing=TRUE),]
                  g <- ggplot(data = df.met, aes(x = submitter, y = metric_value, shape = submission,
                                                 color = submission))
                  submissions <- sort(unique(df.met$submission), decreasing = FALSE)
                  submission.colors <- c("green", "blue", "black", "yellow", "red")[1:length(submissions)]
                  g <- g + scale_color_manual(breaks = submissions,
                                              values = submission.colors)
                  ## g <- g + geom_bar(stat = "identity")
                  g <- g + geom_point(size = 5)
                  g <- g + ylab(paste0("Grand Mean\n(", firstup(met), " Correlation)"))
                  g <- g + xlab(firstup("submitter"))
                  g <- g + coord_flip()
                  g
              })

    gs

}

l_ply(sub.challenges,
      .fun = function(sc) {
          df <- lb.dfs[[sc]]
          gs <- plot.grand.mean(df)
          png(paste0("final-validation-", sc, "-grand-means.png"), width = 480 * 2)
          g <- do.call("grid.arrange", c(gs, nrow = 1, top = paste0(firstup(sc), "-Grained Sub-Challenge")))
          d <- dev.off()
      })

## is_latest is broken. As a temporary work around keep only the best for each submission.
best.lb.dfs <- llply(lb.dfs,
                     .fun = function(df.in) {
                         ## Average over the datasets
                         ## df <- subset(df, ( is_latest == TRUE ) | ( grepl(repo_name, pattern = "baseline")))
                         df.pearson <- subset(df.in, metric == "pearson" & dataset == "Grand mean")
                         df.pearson <- df.pearson[order(df.pearson$metric_value, decreasing=TRUE), ]
                         ##                         df.pearson$tmp_id <- paste0(df.pearson$submitter, df.pearson$repo_name)
                         df.pearson$tmp_id <- paste0(df.pearson$submitter)
                         df.pearson <- df.pearson[!duplicated(df.pearson$tmp_id),]
                         df <- merge(df.in, df.pearson[, c("submitter", "repo_name", "objectId")], all.x = FALSE)
                         print(dim(df))
                         print(unique(df$dataset))
                         df <- subset(df, !grepl(dataset, pattern="mean"))
                         print(dim(df))
                         print(unique(df$dataset))
                         df.ret <- ddply(df, .variables = c("submitter", "repo_name", "objectId", "celltype", "metric"),
                                         .fun = function(df.ct) { data.frame(metric_value = mean(df.ct$metric_value)) })
                         df.ret
                         
                })

latest.lb.dfs <- llply(lb.dfs,
                     .fun = function(df.in) {
                         ## Average over the datasets
                         df <- subset(df.in, submission == "latest")
                         df.ret <- ddply(df, .variables = c("submitter", "repo_name", "objectId", "celltype", "metric"),
                                         .fun = function(df.ct) { data.frame(metric_value = mean(df.ct$metric_value)) })
                         df.ret
                         
                })


plot.cell.type.correlation.heatmap <- function(df, method.id.col = "submitter",
                                               cell.type.col = "celltype",
                                               cor.col = "metric_value",
                                               show.corr.text = FALSE) {

    df <- df[, c(method.id.col, cell.type.col, cor.col)]
    colnames(df) <- c("id", "celltype", "cor")
    df <- df[df$celltype != "Grand mean",]
    df$id <- as.character(df$id)
    df$celltype <- as.character(df$celltype)

    na.rm <- FALSE
    
    cell.type.means <- ddply(df, .variables = c("celltype"),
                             .fun = function(tmp) data.frame(id = "mean", cor = mean(tmp[, "cor"], na.rm=TRUE)))
    cell.type.means <- cell.type.means[order(cell.type.means$cor),]
    
    method.means <- ddply(df, .variables = c("id"),
                   .fun = function(tmp) data.frame(celltype = "mean", cor = mean(tmp[, "cor"], na.rm=na.rm)))
    method.means <- method.means[order(method.means$cor),]

    cell.type.levels <- c(cell.type.means[cell.type.means[, "celltype"] != "mean", "celltype"], "mean")
    id.levels <- c("mean", method.means[method.means[, "id"] != "mean", "id"])

    df <- rbind(df, cell.type.means[, c("id", "celltype", "cor")])
    df <- rbind(df, method.means[, c("id", "celltype", "cor")])    
    df$celltype <- factor(df$celltype, levels = cell.type.levels)
    df$id <- factor(df$id, levels = id.levels)

    df$cor.label <- formatC(df$cor, format="f", digits=2)
    g <- ggplot(data = df, aes(y = id, x = celltype, fill = cor))
    g <- g + geom_tile()
    if(show.corr.text) {
        g <- g + geom_text(aes(label = cor.label))
    }
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   axis.text.y = element_text(size = 5), title = element_text(size = 8))
    g <- g + ylab("Method") + xlab("")
    ## g <- g + scale_fill_continuous("Pearson\ncorrelation", limits = c(-1,1))
    ## g <- g + scale_fill_gradient2("Pearson\ncorrelation", limits = c(-1,1),
    ##                               low = "red", high = "blue", mid = "white", na.value = "black")
    g <- g + scale_fill_gradient2("Correlation", limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
    ## g <- g + theme(text = element_text(size=20))
    g
}


l_ply(sub.challenges,
      .fun = function(sc) {
          ## df <- best.lb.dfs[[sc]]
          df <- latest.lb.dfs[[sc]]
          df.pearson <- subset(df, metric == "pearson")
          df.spearman <- subset(df, metric == "spearman")
          g.pearson <- plot.cell.type.correlation.heatmap(df.pearson)
          g.pearson <- g.pearson + ggtitle("Pearson Correlation")
          g.spearman <- plot.cell.type.correlation.heatmap(df.spearman)
          g.spearman <- g.spearman + ggtitle("Spearman Correlation")          
          png(paste0("final-validation-", sc, "-cell-type-correlations.png"), width = 480 * 2)
          g <- do.call("grid.arrange", c(list(g.pearson, g.spearman),
                                         nrow = 1, top = paste0(firstup(sc), "-Grained Sub-Challenge")))
          d <- dev.off()
      })

l_ply(sub.challenges,
      .fun = function(sc) {
          sz <- 12
          ## df <- best.lb.dfs[[sc]]
          df <- latest.lb.dfs[[sc]]
          df.pearson <- subset(df, metric == "pearson")
          df.spearman <- subset(df, metric == "spearman")
          g.pearson <- plot.cell.type.correlation.heatmap(df.pearson, show.corr.text = TRUE)
          g.pearson <- g.pearson + scale_fill_gradient2(name = "Pearson\nCorrelation",
                                                        limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
          g.pearson <- g.pearson + theme(axis.text.x = element_text(size = sz),
                                         axis.text.y = element_text(size = sz), title = element_text(size = sz))
          
          g.spearman <- plot.cell.type.correlation.heatmap(df.spearman, show.corr.text = TRUE)
          g.spearman <- g.spearman + scale_fill_gradient2(name = "Spearman\nCorrelation",
                                                          limits = c(-1,1), low = "red", high = "blue", mid = "white", na.value = "black")
          g.spearman <- g.spearman + theme(axis.text.x = element_text(size = sz),
                                           axis.text.y = element_text(size = sz), title = element_text(size = sz))
          
          png(paste0("final-validation-pearson-", sc, "-cell-type-correlations.png"), width = 480 * 2)
          print(g.pearson)
          d <- dev.off()
          png(paste0("final-validation-spearman-", sc, "-cell-type-correlations.png"), width = 480 * 2)
          print(g.spearman)
          d <- dev.off()
      })
