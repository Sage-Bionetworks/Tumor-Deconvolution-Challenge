suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(cowplot))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

source("../utils.R")

set.seed(1234)

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

synLogin()

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## Read in the bootstraps-validation-results.rds file
synId <- "syn22951683"
obj <- synGet(synId, downloadFile=TRUE)
results <- readRDS(obj$path)

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


stop("stop")

suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))

calculate.empirical.bayes <-
    function(df, col.id, numerator.id, denominator.id, sample.id.cols, score.col) {
        flag <- df[, col.id] %in% c(numerator.id, denominator.id)
        tmp <- df[flag, ]
        n.tot <- nrow(tmp) / 2
        n.num <-
            sum(unlist(
                dlply(tmp, .variables = sample.id.cols,
                      .fun = function(df.2) {
                          if(nrow(df.2) != 2) {
                              print(df.2)
                              print(c(numerator.id, denominator.id))
                              stop("Was expecting 2 rows\n")
                          }
                          if(any(is.na(df.2[, score.col]))) { stop("Got NA scores\n") }
                          rownames(df.2) <- df.2[, col.id]
                          diff <- df.2[numerator.id, score.col] -
                              df.2[denominator.id, score.col]
                          if(diff > 0) { return(1) }
                          return(0)
                      })))
        n.num / (n.tot - n.num)
    }

## Modified slightly from
## https://stackoverflow.com/questions/53170465/how-to-make-a-base-r-style-boxplot-using-ggplot2
geom_boxplotMod <- function(mapping = NULL, data = NULL, stat = "boxplot", 
    position = "dodge2", ..., outlier.colour = NULL, outlier.color = NULL, 
    outlier.fill = NULL, outlier.shape = 1, outlier.size = 1.5, 
    outlier.stroke = 0.5, outlier.alpha = NULL, notch = FALSE, notchwidth = 0.5,
    varwidth = FALSE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
    linetype = "dashed") # to know how these come here use: args(geom_boxplot)
    {
    list(geom_boxplot(
            mapping = mapping, data = data, stat = stat, position = position,
            outlier.colour = outlier.colour, outlier.color = outlier.color, 
            outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
            outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
            outlier.alpha = outlier.alpha, notch = notch, 
            notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
            show.legend = show.legend, inherit.aes = inherit.aes, linetype = 
            linetype, ...),
        stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.25),
        #the width of the error-bar heads are decreased
        stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.25),
        stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), ...)
        )
    }



        ret.list <- list("bootstrapped.cors" = bootstrapped.cors,
                         "bootstrapped.scores" = bootstrapped.scores,
                         "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "means.by.cell.type.method" = means.by.cell.type.method,
                         "means.over.dataset" = means.over.dataset,
                         "top.performers" = top.performers,
                         "bayes.factors" = bayes.factors)

make.round.text <- function(round, round.str = "Round") {
    round.text <- ""
    if(round == "latest") {
        round.text <- paste0("Latest ", round.str)
    } else if (round == "1") {
        round.text <- paste0(round.str, " 1")
    } else {
        round.text <- paste0("Latest ", round.str, " up to ", round.str, " ", round)
    }
    round.text
}


plot.bootstrap.analysis <-
    function(bootstrapped.scores, mean.bootstrapped.scores,
             means.by.cell.type.method,
             means.over.dataset,
             postfix) {

        cat(paste0("Calculating boxplots\n"))
        boxplots <- list()
        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            mean.scores <- mean.bootstrapped.scores[[sub.challenge]]

            o <- order(mean.scores$pearson)
            mean.scores <- mean.scores[o, ]
            scores[, method.name.col] <- factor(scores[, method.name.col], levels = mean.scores[, method.name.col])
            scores <- na.omit(scores)
            
            g1 <- ggplot(data = scores)
            g1 <- g1 + geom_boxplot(aes_string(x = method.name.col, y = "pearson"))
            g1 <- g1 + coord_flip()
            g1 <- g1 + xlab("Method")
            g1 <- g1 + ylab("Pearson Correlation")
            g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))
            
            g2 <- ggplot(data = scores)
            g2 <- g2 + geom_boxplot(aes_string(x = method.name.col, y = "spearman"))
            g2 <- g2 + coord_flip()
            g2 <- g2 + xlab("Method")
            g2 <- g2 + ylab("Spearman Correlation")
            g2 <- g2 + theme(text = element_text(size=18))    

            boxplots[[paste0(sub.challenge,  "-pearson")]] <- g1
            boxplots[[paste0(sub.challenge,  "-spearman")]] <- g2            

            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)
            title <- paste0(title, " (", round.text, ")")
            ## png(paste0(figs.dir, "rerun-validation-score-boxplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## g <- grid.arrange(g1, g2, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))
            ## d <- dev.off()
        }

        barplots <- list()
        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            mean.scores <- mean.bootstrapped.scores[[sub.challenge]]

            o <- order(mean.scores$pearson)
            mean.scores <- mean.scores[o, ]
            scores[, method.name.col] <- factor(scores[, method.name.col], levels = mean.scores[, method.name.col])
            scores <- na.omit(scores)
            mean.scores[, method.name.col] <- factor(mean.scores[, method.name.col], levels = mean.scores[, method.name.col])
            mean.scores <- na.omit(mean.scores)
            
            g1 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson"))
            g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
            g1 <- g1 + coord_flip()
            g1 <- g1 + xlab("Method")
            g1 <- g1 + ylab("Pearson Correlation")
            g1 <- g1 + ylim(c(-0.25, 1))
            g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))

            g2 <- ggplot(data = mean.scores)
            g2 <- g2 + geom_col(aes_string(x = method.name.col, y = "spearman"), fill = "#E69F00")
            g2 <- g2 + coord_flip()
            g2 <- g2 + xlab("Method")
            g2 <- g2 + ylab("Spearman Correlation")
            g2 <- g2 + ylim(c(-0.25, 1))            
            g2 <- g2 + theme(text = element_text(size=18))    
            g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.ticks.y = element_blank())

            barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
            barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            

            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)            
            title <- paste0(title, " (", round.text, ")")
            ## png(paste0(figs.dir, "rerun-validation-score-box-and-barplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## grid.arrange(g1, g2, nrow=1, widths = c(3, 1), top = textGrob(title, gp = gpar(fontsize = 25)))
            ## d <- dev.off()
        }
        
        ## Spot check that the methods have the same scores for coarse- and fine-grained
        ## NB: some methods may differ between coarse and fine-grained; pick several
        ## baseline/comparator methods that we know to be the same
        check.methods <- c("CIBERSORT", "MCP-counter", "CIBERSORTx")
        check.methods <- sort(unique(c(as.character(means.by.cell.type.method[["coarse"]][["pearson"]][, method.name.col]),
                                       as.character(means.by.cell.type.method[["fine"]][["pearson"]][, method.name.col]))))
        for(meth in check.methods) {

            res.coarse <- means.by.cell.type.method[["coarse"]][["pearson"]]
            res.fine <- means.by.cell.type.method[["fine"]][["pearson"]]            
            
            meth.res.coarse <- res.coarse[res.coarse[, method.name.col] == meth, ]
            meth.res.fine <- res.fine[res.fine[, method.name.col] == meth, ]            
            m <- merge(meth.res.coarse, meth.res.fine, by = c(cell.type.col))
            if(nrow(m) == 0) { next }
            m$diff <- m$cor.x - m$cor.y
            m <- m[!is.na(m$cor.x),]
            eps <- 10^-4
            cat(paste0(meth, " ", postfix, " max diff between coarse and fine-grained is: ", max(abs(m$diff)), "\n"))
            flag <- abs(m$diff) > eps 
            if(any(flag)) {
                print(head(m[flag,,drop=F]))
                cat(paste0("Max diff exceeded for ", meth, " ", postfix, ": ", max(abs(m$diff)), "\n"))
            }
        }

        cat(paste0("Plotting strip plots\n"))
        strip.plots <- list()
        for(sub.challenge in sub.challenges) {
            means <- means.over.dataset[[sub.challenge]][["pearson"]]
            
            g <- plot.strip.plots(means, id.var = method.name.col, cell.type.var = cell.type.col, var = "cor")
            strip.plots[[sub.challenge]] <- g
            
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)
            title <- paste0(title, " (", round.text, ")")
            g <- g + ggtitle(title)

            ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-strip-plot-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## print(g)
            ## d <- dev.off()
        }
        
        cat(paste0("Plotting heatmaps\n"))
        heatmaps <- list()
        for(sub.challenge in sub.challenges) {
            means <- means.by.cell.type.method[[sub.challenge]][["pearson"]]
            
            g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                    id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
##            g <- plot.cell.type.correlation.strip.plots(means, show.corr.text = TRUE, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")

            heatmaps[[sub.challenge]] <- g
            
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)
            title <- paste0(title, " (", round.text, ")")
            
            g <- g + ggtitle(title)
            ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## print(g)
            ## d <- dev.off()
        }

        ret.list <- list("boxplots" = boxplots,
                         "barplots" = barplots,
                         "strip.plots" = strip.plots,
                         "heatmaps" = heatmaps)

        ret.list
        
}

## for(round in c("1", "2", "3", "latest")) {
plots <- list()
for(round in c("2", "1", "3", "latest")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    bootstrapped.cors <- results[[round]][["bootstrapped.cors"]]
    bootstrapped.scores <- results[[round]][["bootstrapped.scores"]]
    mean.bootstrapped.scores <- results[[round]][["mean.bootstrapped.scores"]]
    means.by.cell.type.method <- results[[round]][["means.by.cell.type.method"]]
    means.over.dataset <- results[[round]][["means.over.dataset"]]
    top.performers <- results[[round]][["top.performers"]]
    bayes.factor <- results[[round]][["bayes.factors"]]

    plots[[round]] <- plot.bootstrap.analysis(bootstrapped.scores, mean.bootstrapped.scores,
                                              means.by.cell.type.method,
                                              means.over.dataset,
                                              postfix)
    
    cat("Bayes factor resultes K <= 3 (or K <= 5) suggests a tie\n")
    for(sub.challenge in sub.challenges) {
        best.team <- results[[round]][["top.performers"]][[sub.challenge]]
        cat(paste0("Top performer for Round ", round, " in ", sub.challenge, " sub-challenge: ",
                   best.team, "\n"))
        tbl <- results[[round]][["bayes.factors"]][[sub.challenge]]
        colnames(tbl)[1] <- "team"
        print(tbl)
        file <- paste0("rerun-validation-bayes-factors-round-", round, "-sc-", sub.challenge, "-vs-",
                       make.names(best.team), ".tsv")
        write.table(file = file, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    }
}

plot.scores.over.rounds <- function(df) {
    df.latest <- subset(df, Round == "latest")
    o <- order(df.latest$pearson)
    lvls <- df.latest[o, method.name.col]
    df <- subset(df, Round != "latest")
    df[, method.name.col] <- factor(df[, method.name.col], levels = lvls)
    df$Submission <- df$Round
    g <- ggplot()
    g <- g + geom_point(data = df, aes_string(x = "pearson", y = method.name.col, colour = "Submission"))
    g <- g + xlab("Pearson Correlation") + ylab("Method")
    g
}

all.means <-
    llply(sub.challenges,
          .fun = function(subchallenge) {
              ret <- ldply(results, .fun = function(df) df$mean.bootstrapped.scores[[subchallenge]])
              colnames(ret)[1] <- "Round"
              ret <- na.omit(ret)
              ## Round x score means score at round x or the latest round < x.
              ## i.e., there may be duplicate results in this table.
              ## Rather than using duplicated(pearson) to remove them, just
              ## order by round and they will be overlaid on the plot below.
              o <- order(ret$Round, decreasing = TRUE)
              ret <- ret[o,]
              ret
          })

g.score.vs.round <-
    llply(sub.challenges, .fun = function(subchallenge) plot.scores.over.rounds(all.means[[subchallenge]]))

if(FALSE) {
    l_ply(sub.challenges,
          .fun = function(subchallenge) {
              g <- g.score.vs.round[[subchallenge]]
              png(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", subchallenge, ".png"))
              print(g)
              d <- dev.off()
          })
}

            barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
            barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            


g.bootstrap.coarse.pearson.round1 <- plots[["1"]][["barplots"]][["coarse-pearson"]]
g.bootstrap.coarse.spearman.round1 <- plots[["1"]][["barplots"]][["coarse-spearman"]]

title <- "Coarse-Grained Sub-Challenge (First Submission)"
g.bootstrap.coarse.round1 <- grid.arrange(g.bootstrap.coarse.pearson.round1,
                                          g.bootstrap.coarse.spearman.round1, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))

g.bootstrap.fine.pearson.round1 <- plots[["1"]][["barplots"]][["fine-pearson"]]
g.bootstrap.fine.spearman.round1 <- plots[["1"]][["barplots"]][["fine-spearman"]]

title <- "Fine-Grained Sub-Challenge (First Submission)"
g.bootstrap.fine.round1 <- grid.arrange(g.bootstrap.fine.pearson.round1,
                                        g.bootstrap.fine.spearman.round1, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))

g.round.coarse <- g.score.vs.round[["coarse"]]
g.round.coarse <- g.round.coarse + ggtitle("Coarse-Grained Sub-Challenge")

g.round.fine <- g.score.vs.round[["fine"]]
g.round.fine <- g.round.fine + ggtitle("Fine-Grained Sub-Challenge")

g <- plot_grid(g.bootstrap.coarse.round1, g.bootstrap.fine.round1, g.round.coarse, g.round.fine)

png(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", "all", postfix, ".png"), width = 2 * 480)                    
grid.arrange(g1, g2, nrow=1)
d <- dev.off()

  - bootstraps for coarse (A) and fine (B) (round 1)
  - performance over rounds for coarse (C) and fine (D)


cat("Exiting successfully\n")
q(status=0)
