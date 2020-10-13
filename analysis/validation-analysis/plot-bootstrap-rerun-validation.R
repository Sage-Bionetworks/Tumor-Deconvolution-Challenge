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
    function(res, bootstrapped.scores, mean.bootstrapped.scores,
             means.by.cell.type.method,
             means.over.dataset,
             postfix) {

        top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
        priority.methods <- unique(c(top.performers, unique(subset(res, comparator==TRUE)[, method.name.col])))

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
            
            ## title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            ## round.text <- make.round.text(round)
            ## title <- paste0(title, " (", round.text, ")")
            ## g <- g + ggtitle(title)

            ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-strip-plot-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## print(g)
            ## d <- dev.off()

            flag <- means[, method.name.col] %in% priority.methods
            g <- plot.strip.plots(means[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor")
            strip.plots[[paste0(sub.challenge, "-priority")]] <- g
            
        }

        coarse.means <- means.over.dataset[["coarse"]][["pearson"]]
        fine.means <- means.over.dataset[["fine"]][["pearson"]]        
        all.means <- rbind(coarse.means, fine.means)
        
        all.means <-
            ddply(all.means, .variables = c(method.name.col, cell.type.col, "boot.i"),
                  .fun = function(df) {
                      data.frame(cor = mean(df$cor))
                  })
        
        
        g <- plot.strip.plots(all.means, id.var = method.name.col, cell.type.var = cell.type.col, var = "cor")
        strip.plots[["merged"]] <- g

        flag <- all.means[, method.name.col] %in% priority.methods
        g <- plot.strip.plots(all.means[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor")
        strip.plots[["merged-priority"]] <- g

        
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
        
        coarse.means <- means.by.cell.type.method[["coarse"]][["pearson"]]
        fine.means <- means.by.cell.type.method[["fine"]][["pearson"]]        
        all.means <- rbind(coarse.means, fine.means)
        
        all.means <-
            ddply(all.means, .variables = c(method.name.col, cell.type.col),
                  .fun = function(df) {
                      data.frame(cor = mean(df$cor))
                  })
        g <- plot.cell.type.correlation.heatmap(all.means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        heatmaps[["merged"]] <- g
        
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

    res.round <- results[[round]][["res.round"]]
    bootstrapped.cors <- results[[round]][["bootstrapped.cors"]]
    bootstrapped.scores <- results[[round]][["bootstrapped.scores"]]
    mean.bootstrapped.scores <- results[[round]][["mean.bootstrapped.scores"]]
    means.by.cell.type.method <- results[[round]][["means.by.cell.type.method"]]
    means.over.dataset <- results[[round]][["means.over.dataset"]]
    top.performers <- results[[round]][["top.performers"]]
    bayes.factors <- results[[round]][["bayes.factors"]]

    if(FALSE) {
        ranges <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      tmp <- res.round[res.round[, subchallenge.col] == sub.challenge, ]
                      ddply(tmp, .variables = c(cell.type.col, "distribution.type"),
                            .fun = function(df) {
                                mn <- min(df[, measured.col], na.rm=TRUE)
                                mx <- max(df[, measured.col], na.rm=TRUE)
                                data.frame(min = mn, max = mx, range = mx - mn)
                            })
                  })
    }
                  
    plots[[round]] <- plot.bootstrap.analysis(res.round, bootstrapped.scores, mean.bootstrapped.scores,
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
    order.round <- "1"
    df.order <- subset(df, Round == order.round)
    o <- order(df.order$pearson)
    lvls <- df.order[o, method.name.col]
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

shift.limit <- function(val) {
    sign(val) * round(abs(val) + 0.05, digits = 1)
}

source("make-validation-performance-figs.R")

cat("Exiting successfully\n")
q(status=0)
