suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))

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

## Read in the rerun predicitons (i.e., where the coarse- and fine-grained datasets are the same)
synId <- "syn22320329"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

## Read in the bootstraps.rds file
synId <- "syn22344963"
obj <- synGet(synId, downloadFile=TRUE)
bootstraps <- readRDS(obj$path)

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
tmp <- unique(res.all[, !(colnames(res.all) %in% c(cell.type.col, prediction.col, measured.col))])
cell.types.by.sub.challenge <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col)])
tmp <- merge(tmp, cell.types.by.sub.challenge, all = TRUE)
measured.tbl <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, measured.col)])
prediction.tbl <- res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, method.name.col, round.col, prediction.col)]
tmp <- merge(tmp, measured.tbl, all = TRUE)
tmp <- merge(tmp, prediction.tbl, all = TRUE)
res.all <- tmp
flag <- !is.na(res.all[, measured.col])
res.all <- res.all[flag, ]

res.all[, cell.type.col] <- as.character(res.all[, cell.type.col])
res.all <- rename.cell.types(res.all, from.col = cell.type.col, to.col = cell.type.col)

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

do.bootstrap.analysis <-
    function(res.input, bootstraps, 
             method.name.col, 
             model.id.col, subchallenge.col, measured.col, cell.type.col,
             dataset.name.col, sample.id.col, prediction.col,
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
        
        res.round <- merge(res.input, submitter.tbl, by = c(method.name.col, subchallenge.col, round.col))
        
        method.id.col <- method.name.col
        
        res <- res.round
        
        tbls <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      flag <- res[, subchallenge.col] == sub.challenge
                      tmp <- res[flag, ]
                      tmp
                  })
        
        tbls.by.cell <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      tmp <- tbls[[sub.challenge]]
                      flag <- !is.na(tmp[, measured.col])
                      tmp <- tmp[flag, ]
                      cells <- unique(tmp[, cell.type.col])
                      names(cells) <- cells
                      llply(cells,
                            .fun = function(ct) {
                                flag <- tmp[, cell.type.col] == ct
                                ret <- tmp[flag, ]
                                ret$id <- paste0(ret[, dataset.name.col], "-", ret[, sample.id.col])
                                ret
                            })
                  })

        ## Calculate both pearson and spearman correlation over bootstraps, within dataset and cell type
        cat(paste0("Calculating bootstrap correlations\n"))
        bootstrapped.cors <-
            llply(sub.challenges, 
                  .fun = function(sub.challenge) {
                      tmp <- tbls[[sub.challenge]]                      
                      flag <- !is.na(tmp[, measured.col])
                      tmp <- tmp[flag, ]
                      methods <- unique(tmp[, method.id.col])
##                      methods <- c("Aginome-XMU", "CIBERSORTx")
                      names(methods) <- methods
                      n.bootstraps <- length(bootstraps[[sub.challenge]])
                      indices <- 1:n.bootstraps
                      names(indices) <- indices
                      ret.all <-
                          ldply(methods, .parallel = TRUE,
                                .fun = function(method) {
                                    print(method)
                                    ret.i <-
                                        ldply(indices,
                                              .fun = function(i) {
                                                  ret.method <-
                                                      ldply(tbls.by.cell[[sub.challenge]],
                                                            .fun = function(df.in) {
                                                                flag <- df.in[, method.id.col] == method
                                                                df <- df.in[flag, ]
                                                                if(any(duplicated(df$id))) {
                                                                    print(head(df[my.dup(df$id),]))
                                                                    stop("stop")
                                                                }
                                                                rownames(df) <- df$id
                                                                sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                                                sample.ids <- sample.ids[sample.ids %in% df$id]
                                                                if(!(all(sample.ids %in% rownames(df)))) {
                                                                    stop("Some sample ids not in df\n")
                                                                }
                                                                df <- df[sample.ids,]
                                                                df
                                                            })
                                                  ret.method <- ret.method[, -1]
                                                  score <-
                                                      ddply(ret.method, .variables = c(dataset.name.col),
                                                            .fun = function(df.ds) {
                                                                tmp <- ddply(df.ds, .variables = c(cell.type.col),
                                                                             .fun = function(df.ct) {
                                                                                 if(any(is.na(df.ct[,2]))) {
                                                                                     print(df.ct)
                                                                                 }
                                                                                 mts <- c("pearson", "spearman", "rmse")
                                                                                 names(mts) <- mts
                                                                                 vals <- llply(mts,
                                                                                               .fun = function(comp.method) {
                                                                                                   pred <- as.numeric(df.ct[, prediction.col])
                                                                                                   measured <- as.numeric(df.ct[, measured.col])
                                                                                                   
                                                                                                   if(any(is.na(pred))) {
                                                                                                       if(!all(is.na(pred))) {
                                                                                                           stop("Some but not all NA\n")
                                                                                                       }
                                                                                                       return(NA)
                                                                                                   }
                                                                                                   if(var(pred) == 0) { return(0) }

                                                                                                   val <- NA
                                                                                                   if(comp.method %in% c("pearson", "spearman")) {
                                                                                                       val <- cor(pred, measured,
                                                                                                                  method = comp.method)
                                                                                                   } else {
                                                                                                       val <- sqrt(mean((pred-measured)^2))
                                                                                                   }
                                                                                                   val
                                                                                               })
                                                                                 as.data.frame(vals)
                                                                             })
                                                                colnames(tmp)[1] <- cell.type.col
                                                                tmp
                                                            })
                                                  score
                                              })
                                    colnames(ret.i)[1] <- "boot.i"
                                    ret.i
                                })
                      colnames(ret.all)[1] <- method.id.col
                      ret.all
                  })
        print(colnames(bootstrapped.cors[[1]]))
        print(head(bootstrapped.cors[[1]]))

        cat(paste0("Calculating bootstrapped scores\n"))
        bootstrapped.scores <-
            llply(sub.challenges, .parallel = TRUE,
                  .fun = function(sub.challenge) {
                      df <- bootstrapped.cors[[sub.challenge]]
                      flag <- res.round[, subchallenge.col] == sub.challenge
                      un <- unique(res.round[flag, c(method.id.col, method.name.col)])
                      df <- merge(df, un)
                      ## Average over cell type (within method and dataset and bootstrap sample)
                      df <- ddply(df,
                                  .variables = c(method.id.col, dataset.name.col, "boot.i"),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                                  })
                      ## Average over dataset (within method and bootstrap sample)
                      df <- ddply(df,
                                  .variables = c(method.id.col, "boot.i"),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                                  })
                      df
                  })
    
        cat(paste0("Calculating mean bootstrapped scores\n"))
        mean.bootstrapped.scores <-
            llply(bootstrapped.scores, .parallel = TRUE,
                  .fun = function(df) {
                      ## Average over bootstraps
                      df <- ddply(df,
                                  .variables = c(method.id.col),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                                  })
                      o <- order(df$pearson, decreasing = TRUE)
                      df <- df[o,]
                      df
                  })

        cat(paste0("Calculating mean by cell type\n"))        
        means.by.cell.type.method <-
            llply(bootstrapped.cors,
                  .fun = function(df) {
                      methods <- c("pearson", "spearman", "rmse")
                      na.rm <- FALSE
                      names(methods) <- methods
                      res <- llply(methods,
                                   .fun = function(method) {
                                       ## first, average over bootstrap
                                       ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                                    .fun = function(df) {
                                                        data.frame(cor = mean(df[, method], na.rm=na.rm))
                                                    })
                                       ## now, average over dataset
                                       ret2 <- ddply(ret, .variables = c(method.id.col, cell.type.col),
                                                     .fun = function(df) {
                                                         data.frame(cor = mean(df$cor, na.rm=na.rm))
                                                     })
                                   })
                  })

        means.over.dataset <-
            llply(bootstrapped.cors,
                  .fun = function(df) {
                      methods <- c("pearson", "spearman", "rmse")
                      na.rm <- FALSE
                      names(methods) <- methods
                      res <- llply(methods,
                                   .fun = function(method) {
                                       ## average over dataset
                                       ret <- ddply(df, .variables = c(method.id.col, cell.type.col, "boot.i"),
                                                    .fun = function(df) {
                                                        data.frame(cor = mean(df[, method], na.rm=na.rm))
                                                    })
                                   })
                  })

        cat(paste0("Calculating bayes factors\n"))
        top.performers <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
                      mean.scores <- na.omit(mean.scores)
                      scores <- bootstrapped.scores[[sub.challenge]]
                      numerator.indx <- 1
                      numerator.id <- mean.scores[numerator.indx, method.id.col]
                      numerator.id
                  })
        
        bayes.factors <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
                      mean.scores <- na.omit(mean.scores)
                      scores <- bootstrapped.scores[[sub.challenge]]
                      numerator.indx <- 1
                      numerator.id <- mean.scores[numerator.indx, method.id.col]
                      indices <- 1:nrow(mean.scores)
                      indices <- indices[!(indices == numerator.indx)]
                      names(indices) <- mean.scores[indices, method.id.col]
                      ret <-
                          ldply(indices,
                                .fun = function(i) {
                                    denominator.id <- mean.scores[i, method.id.col]
                                    methods <- c("pearson", "spearman")
                                    names(methods) <- methods
                                    res <- llply(methods,
                                                 .fun = function(method) {
                                                     bf <- calculate.empirical.bayes(scores, col.id = method.id.col,
                                                                                     numerator.id = numerator.id,
                                                                                     denominator.id = denominator.id,
                                                                                     sample.id.cols = "boot.i",
                                                                                     score.col = method)
                                                 })
                                    ret.df <- as.data.frame(res)
                                    ret.df
                                })
                      ret
                  })

        ret.list <- list("bootstrapped.cors" = bootstrapped.cors,
                         "bootstrapped.scores" = bootstrapped.scores,
                         "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "means.by.cell.type.method" = means.by.cell.type.method,
                         "means.over.dataset" = means.over.dataset,
                         "top.performers" = top.performers,
                         "bayes.factors" = bayes.factors)
                         
        return(ret.list)

        
        cat(paste0("Calculating boxplots\n"))
        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, method.name.col, drop=FALSE])
            mean.scores <- merge(mean.scores, un)
            scores <- merge(scores, un)    
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

            png(paste0(figs.dir, "rerun-validation-score-boxplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- ""
            if(round == "latest") {
                round.text <- "Latest Round"
            } else if (round == "1") {
                round.text <- "Round 1"
            } else {
                round.text <- paste0("Latest Round up to Round ", round)
            }
            title <- paste0(title, " (", round.text, ")")
            grid.arrange(g1, g2, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))
            d <- dev.off()
        }

        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, method.name.col, drop=FALSE])
            mean.scores <- merge(mean.scores, un)
            scores <- merge(scores, un)    
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

            png(paste0(figs.dir, "rerun-validation-score-box-and-barplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- ""
            if(round == "latest") {
                round.text <- "Latest Round"
            } else if (round == "1") {
                round.text <- "Round 1"
            } else {
                round.text <- paste0("Latest Round up to Round ", round)
            }
            title <- paste0(title, " (", round.text, ")")
            grid.arrange(g1, g2, nrow=1, widths = c(3, 1), top = textGrob(title, gp = gpar(fontsize = 25)))
            d <- dev.off()
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
        for(sub.challenge in sub.challenges) {
            means <- means.over.dataset[[sub.challenge]][["pearson"]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, method.name.col, drop=FALSE])
            means <- merge(means, un)
            
            g <- plot.strip.plots(means, id.var = method.name.col, cell.type.var = cell.type.col, var = "cor")
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- ""
            if(round == "latest") {
                round.text <- "Latest Round"
            } else if (round == "1") {
                round.text <- "Round 1"
            } else {
                round.text <- paste0("Latest Round up to Round ", round)
            }
            title <- paste0(title, " (", round.text, ")")
            
            g <- g + ggtitle(title)
            png(paste0(figs.dir, "rerun-validation-bootstrap-cell-strip-plot-", sub.challenge, postfix, ".png"), width = 2 * 480)
            print(g)
            d <- dev.off()

        }
        
        cat(paste0("Plotting heatmaps\n"))
        for(sub.challenge in sub.challenges) {
            means <- means.by.cell.type.method[[sub.challenge]][["pearson"]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, method.name.col, drop=FALSE])
            means <- merge(means, un)
            
            g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                    id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
##            g <- plot.cell.type.correlation.strip.plots(means, show.corr.text = TRUE, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- ""
            if(round == "latest") {
                round.text <- "Latest Round"
            } else if (round == "1") {
                round.text <- "Round 1"
            } else {
                round.text <- paste0("Latest Round up to Round ", round)
            }
            title <- paste0(title, " (", round.text, ")")
            
            g <- g + ggtitle(title)
            png(paste0(figs.dir, "rerun-validation-bootstrap-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
            print(g)
            d <- dev.off()
        }

        ret.list <- list("bootstrapped.cors" = bootstrapped.cors,
                         "bootstrapped.scores" = bootstrapped.scores,
                         "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "top.performers" = top.performers,
                         "bayes.factors" = bayes.factors,
                         "means.by.cell.type.method" = means.by.cell.type.method)
        ret.list
        
}

results <- list()
## for(round in c("1", "2", "3", "latest")) {
for(round in c("2", "1", "3", "latest")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    if(FALSE) {
        res.input <- res.all
        round.col <- "submission"
        round <- "2"
        postfix <- paste0("-round-", round)
    }
    
    results[[round]] <- do.bootstrap.analysis(res.all, bootstraps, 
                                              method.name.col,
                                              model.id.col, subchallenge.col, measured.col, cell.type.col,
                                              dataset.name.col, sample.id.col, prediction.col,
                                              round.col = "submission", round = round,
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

save.image(".Rdata.plot.bootstrap")

## Store the above results to Synapse
rounds <- names(results)
    
parent.id <- "syn22320184"

url.base <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
this.script <- "perform-bootstrap-rerun-validation-analysis.R"
script_url <- paste0(url.base, "/", "validation-analysis", "/", this.script)

file <- "bootstrap-validation-results.rds"
saveRDS(results, file)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = parent.id, executed = script_url, synapseStore = TRUE)
synStore(f)

if(FALSE) {
    for(round in rounds) {
        res <- results[[round]]
        for(nm in names(res)) {
            file <- paste0(nm, ".tsv")
            write.table(file = file, res[[nm]], col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
            cat(paste0("Storing ", file, " to synapse\n"))
            f <- File(file, parentId = parent.id, executed = script_url, synapseStore = TRUE)
            synStore(f)
        }
    }
}

subchallenges <- c("coarse", "fine")
names(subchallenges) <- subchallenges

plot.scores.over.rounds <- function(df) {
    df.latest <- subset(df, Round == "latest")
    o <- order(df.latest$pearson)
    lvls <- df.latest[o, "method.name"]
    df <- subset(df, Round != "latest")
    df$method.name <- factor(df$method.name, levels = lvls)
    g <- ggplot()
    g <- g + geom_point(data = df, aes(x = pearson, y = method.name, colour = Round))
    g <- g + xlab("Pearson Correlation") + ylab("Method")
    g
}

all.means <-
    llply(subchallenges,
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
    llply(subchallenges, .fun = function(subchallenge) plot.scores.over.rounds(all.means[[subchallenge]]))

l_ply(subchallenges,
      .fun = function(subchallenge) {
          g <- g.score.vs.round[[subchallenge]]
          png(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", subchallenge, ".png"))
          print(g)
          d <- dev.off()
      })

g1 <- g.score.vs.round[["coarse"]]
g1 <- g1 + ggtitle("Coarse-Grained Sub-Challenge")
g2 <- g.score.vs.round[["fine"]]
g2 <- g2 + ggtitle("Fine-Grained Sub-Challenge")

png(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", "all", postfix, ".png"), width = 2 * 480)                    
grid.arrange(g1, g2, nrow=1)
d <- dev.off()

save.image(".Rdata.plot.bootstrap")

stop("stop")

cat("Exiting successfully\n")
q(status=0)