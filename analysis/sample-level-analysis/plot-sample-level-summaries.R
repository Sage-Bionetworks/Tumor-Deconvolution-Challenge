suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(ggbeeswarm))
suppressPackageStartupMessages(p_load(cowplot)) # for plot_grid

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

source("../utils.R")

set.seed(1234)

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

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"

dataset.annotation <- unique(res.all[, c(dataset.name.col, "tumor.type", "distribution.type", "mixture.type")])

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

## Restrict to "deconvolution" methods -- i.e., those that report fractions or scores that
## can be compared across cell types. Note that in some cases, e.g., Aginome-XMU, one submission
## allows this while another does not. Aginome-XMU uses an ensemble that includes MCP in submissions
## after the first, which precludes its use here. "NA" indicates the method may be used to
## compare across cell types in any round; otherwise the particular round of applicability
## is provided
## RGEGEN_LAB is so bad it is skewing the results; leave it off.
deconv.fraction.methods <-
    list(
        ## "REGGEN_LAB" = NA,
         "NPU" = NA,
         "CCB" = NA,
         "DA_505" = NA,
         "Aginome-XMU" = "1",
         "D3Team" = NA,
         "IZI" = NA,
         "LeiliLab" = NA,         
         "CIBERSORT" = NA,
         "CIBERSORTx" = NA,
         "quanTIseq" = NA,
         "EPIC" = NA,
         "TIMER" = NA)
deconv.score.methods <-
    list("NYIT_glomerular" = NA, "Biogem" = NA, "Patrick" = NA, "TJU" = NA)

deconv.summary <-
    rbind(data.frame(method = names(deconv.fraction.methods), submission = unlist(deconv.fraction.methods), output = "fraction",
                     stringsAsFactors = FALSE),
          data.frame(method = names(deconv.score.methods), submission = unlist(deconv.score.methods), output = "score",
                     stringsAsFactors = FALSE))
          

flags <-
    llply(1:nrow(deconv.summary),
          .fun = function(i) {
              flag <- res.all[, method.name.col] == deconv.summary[i, "method"]
              if(!is.na(deconv.summary[i, "submission"])) {
                  flag <- flag & ( res.all[, round.col] == deconv.summary[i, "submission"])
              }
              flag
          })
flag <- Reduce("|", flags)

res.all <- res.all[flag,]

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")
comparison.metrics <- list("Pearson" = "pearson", "Spearman" = "spearman", "RMSE" = "RMSE")

## A method must report all cell types in a given sub-challenge to be considered in this cross cell-type comparison
reports.all.cell.types <-
    ldply(sub.challenges,
          .fun = function(sub.challenge) {
              flag <- res.all[, subchallenge.col] == sub.challenge
              res.sc <- res.all[flag, ]
              ddply(res.sc, .variables = c(method.name.col),
                    .fun = function(df) {
                        data.frame(reports.all = !any(is.na(df$prediction)))
                    })
          })
colnames(reports.all.cell.types)[1] <- subchallenge.col
reports.all.cell.types <- subset(reports.all.cell.types, reports.all == TRUE)

res.all <- merge(res.all, reports.all.cell.types)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(reshape2))

do.sample.level.analysis <-
    function(res.input, dataset.annotation,
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
        
        ## Calculate metrics across cell types (within a dataset)
        cat(paste0("Calculating metrics\n"))
        metric.all.res <-
            llply(tbls, 
                  .fun = function(tbl) {
                      all.res <-
                          ddply(tbl, .variables = c(method.name.col, dataset.name.col, sample.id.col),
                                .fun = function(df) {
                                    res <- llply(comparison.metrics,
                                                 .fun = function(metric) {
                                                     val <- NA
                                                     pred <- as.numeric(df[, prediction.col])
                                                     obs <- as.numeric(df[, measured.col])
                                                     if(metric == "RMSE") {
                                                         ## Only calculate RMSE for the fraction-based approaches
                                                         if(df[1, method.name.col] %in% names(deconv.fraction.methods)) {
                                                             val <- sqrt(mean((pred - obs)^2))
                                                         } 
                                                     } else {
                                                         val <- cor(pred, obs, method = metric)
                                                     }
                                                     val
                                                 })
                                    ret <- unlist(res)
                                    names(ret) <- firstup(names(res))
                                    ret
                                })
                      all.res
                  })

        ## Average over the samples within a dataset, then over datasets
        metric.sum.res <-
            llply(metric.all.res,
                  .fun = function(tbl) {
                      ## Average over samples within a dataset
                      ret.ds <-
                          ddply(tbl, .variables = c(method.name.col, dataset.name.col), 
                                .fun = function(df) {
                                    res <- llply(comparison.metrics,
                                                 .fun = function(metric) {
                                                     val <- mean(df[, firstup(metric)])
                                                 })
                                    ret <- unlist(res)
                                    names(ret) <- firstup(names(res))
                                    ret
                                })
                      ## Average over datasets
                      ddply(ret.ds, .variables = c(method.name.col),
                            .fun = function(df) {
                                res <- llply(comparison.metrics,
                                             .fun = function(metric) {
                                                 val <- mean(df[, firstup(metric)])
                                             })
                                ret <- unlist(res)
                                names(ret) <- firstup(names(res))
                                ret
                            })
                      })

        g.summaries <-
            llply(metric.sum.res,
                  .fun = function(tbl) {
                      m.tbl <- melt(tbl)
                      o <- order(tbl$Pearson, decreasing = TRUE)
                      lvls <- tbl[o, method.name.col]
                      colnames(m.tbl) <- c(method.name.col, "metric", "val")
                      m.tbl[, method.name.col] <- factor(m.tbl[, method.name.col], levels = lvls)
                      g <- ggplot()
                      g <- g + geom_point(data = m.tbl, aes_string(x = method.name.col, y = "val"))
                      g <- g + facet_wrap("metric", scales = "free_y")
                      g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                      g <- g + ylab("") + xlab("")
                      g 
                  })

        ## Plot beeswarm of individuals, color coded by dataset
        g.swarm <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                      llply(comparison.metrics,
                            .fun = function(metric) {
                                tbl <- metric.all.res[[sub.challenge]]
                                m.tbl <- melt(tbl)
                                colnames(m.tbl) <- c(method.name.col, dataset.name.col, sample.id.col, "metric", "val")
                                m.tbl <- m.tbl[m.tbl[, "metric"] == firstup(metric), ]
                                m.tbl <- na.omit(m.tbl)
                                m.tbl <- merge(m.tbl, dataset.annotation)
                                decreasing <- TRUE
                                if(metric == "RMSE") { decreasing <- FALSE }
                                o <- order(metric.sum.res[[sub.challenge]][, firstup(metric)], decreasing = decreasing)
                                lvls <- metric.sum.res[[sub.challenge]][o, method.name.col]
                                m.tbl[, method.name.col] <- factor(m.tbl[, method.name.col], levels = lvls)
                                g <- ggplot()
                                g <- g + geom_boxplot(data = m.tbl, aes_string(x = dataset.name.col, y = "val",
                                                                               colour = "mixture.type",
                                                                               linetype = "distribution.type"))
                                g <- g + facet_wrap(method.name.col, scales = "free_y")
                                g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                                g <- g + ylab("") + xlab("")
                                g 
                            })
                  })

        ## Plot the correlations
        plot.correlations <- FALSE
        if(plot.correlations) {
            l_ply(sub.challenges,
                  .fun = function(sub.challenge) {
                      tbl <- tbls[[sub.challenge]]
                      d_ply(tbl, .variables = c(method.name.col),
                            .fun = function(df.meth) {
                                file <- paste0("correlations-", sub.challenge, "-", make.names(df.meth[1, method.name.col]), postfix, ".pdf")
                                pdf(file, onefile = TRUE)
                                d_ply(df.meth,
                                      .variables = c(dataset.name.col),
                                      .fun = function(df.ds) {
                                          g <- ggplot(data = df.ds, aes_string(x = measured.col, y = prediction.col))
                                          g <- g + geom_point()
                                          g <- g + facet_wrap(sample.id.col)
                                          g <- g + geom_smooth(method='lm', formula = y ~ x)
                                          print(g)
                                      })
                                d <- dev.off()
                            })
                  })
        }
        
        ret.list <- list("metric.all.res" = metric.all.res, "metric.sum.res" = metric.sum.res,
                         "g.summaries" = g.summaries, "g.swarm" = g.swarm)
        return(ret.list)

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
    
    results[[round]] <- do.sample.level.analysis(res.all, dataset.annotation,
                                                 method.name.col,
                                                 model.id.col, subchallenge.col, measured.col, cell.type.col,
                                                 dataset.name.col, sample.id.col, prediction.col,
                                                 round.col = "submission", round = round,
                                                 postfix)

    g.sum.coarse <- results[[round]][["g.summaries"]][["coarse"]]
    g.sum.coarse <- g.sum.coarse + ggtitle("Coarse-Grained Sub-Challenge")
    g.sum.fine <- results[[round]][["g.summaries"]][["fine"]]
    g.sum.fine <- g.sum.fine + ggtitle("Fine-Grained Sub-Challenge")
    png(paste0("sample-level-metric-summary", postfix, ".png"))
    g <- plot_grid(g.sum.coarse, g.sum.fine, nrow = 2, labels = "AUTO") 
    print(g)
    ## grid.arrange(g.sum.coarse, g.sum.fine)
    d <- dev.off()
}



cat("Exiting successfully\n")

