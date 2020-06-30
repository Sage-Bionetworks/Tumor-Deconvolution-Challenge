suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))

set.seed(1234)

synLogin()

assign.baseline.names <- function(df, repo.name.col = "repo_name") {
    trans <-
        list("baseline_method1" = "CIBERSORT",
             "baseline_method2" = "MCP-counter",
             "baseline_method3" = "quanTIseq",
             "baseline_method4" = "xCell",
             "baseline_method5" = "EPIC",
             "baseline_method6" = "TIMER",
             "baseline_method7" = "CIBERSORTx")             
    for(nm in names(trans)) {
        flag <- (grepl(df[, repo.name.col], pattern = nm))
        df[flag, repo.name.col] <- trans[[nm]]
    }
    df
}

## Read in the final predictions
synId <- "syn22149603"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
model.id.col <- "modelId"
method.name.col <- "repo_name"

res.all$modelId <- paste0(as.character(res.all$repo_name), "-", as.character(res.all$submitterId), "-", as.character(res.all$objectId))

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
tmp <- unique(res.all[, !(colnames(res.all) %in% c(cell.type.col, prediction.col, measured.col))])
cell.types.by.sub.challenge <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col)])
tmp <- merge(tmp, cell.types.by.sub.challenge, all = TRUE)
measured.tbl <- unique(res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, measured.col)])
prediction.tbl <- res.all[, c(cell.type.col, subchallenge.col, sample.id.col, dataset.name.col, model.id.col, prediction.col)]
tmp <- merge(tmp, measured.tbl, all = TRUE)
tmp <- merge(tmp, prediction.tbl, all = TRUE)
res.all <- tmp

## Check that is.latest is set correctly
## NB: if objectId X < objectId Y then X was submitted before Y
baseline.method.flag <- grepl(res.all$repo_name, pattern="baseline")
## Confirm that each baseline method was only submitted once (hence, we should keep all
## baseline method results)
n.unique.baselines <- nrow(unique(res.all[baseline.method.flag, c("objectId", "repo_name", "subchallenge")]))
if(n.unique.baselines != 2 * length(unique(res.all[baseline.method.flag, "repo_name"]))) {
  stop("Did not see each baseline method submitted exactly twice as we had expected to see\n")
}

latest.objectIds <-
        ddply(unique(res.all[!baseline.method.flag, c("objectId", "subchallenge", "submitterId")]),
              .variables = c("subchallenge", "submitterId"),
              .fun = function(df) data.frame(objectId = max(df$objectId)))

## Ensure that is_latest is properly set.
## i.e., it is TRUE if and only if the corresponding object is in latest.objectIds or
## it is a baseline method
is_latest.objectIds <- unique(subset(res.all[!baseline.method.flag, ], is_latest == TRUE)[, "objectId"])
if(!(all(sort(latest.objectIds$objectId) == sort(is_latest.objectIds)))) {
   stop("is_latest flag is not set consistenly with our understanding of how objectIds are defined\n")
}

res.latest <- res.all[baseline.method.flag | (res.all$is_latest == TRUE), ]

res.latest <- assign.baseline.names(res.latest)
res.latest$modelId <- paste0(as.character(res.latest$repo_name), "-", as.character(res.latest$submitterId), "-", as.character(res.latest$objectId))

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
un <- unique(res.latest[, c("repo_name", "submitterId")])
flag <- my.dup(un$repo_name)
dup.repo.names <- un[flag, "repo_name"]
flag <- res.latest$repo_name %in% dup.repo.names
res.latest[flag, "repo_name"] <- paste0(as.character(res.latest[flag, "repo_name"]), "-", as.character(res.latest[flag, "objectId"]))

res <- res.latest

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

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
              flag <- res[, subchallenge.col] == sub.challenge
              tmp <- res[flag, ]
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
    

n.bootstraps <- 100

bootstraps <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              flag <- res[, subchallenge.col] == sub.challenge
              tmp <- res[flag, ]
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              datasets <-
                  dlply(tmp, .variables = c(dataset.name.col),
                        .fun = function(df) {
                            unique(df[, sample.id.col])
                        })
              dataset.names <- names(datasets)
              names(dataset.names) <- dataset.names
              llply(1:n.bootstraps,
                    .fun = function(i) {
                        ret <- ldply(datasets,
                                     .fun = function(ds) {
                                         tmp.r <- data.frame(sample.id = sample(ds, size = length(ds), replace = TRUE))
                                         colnames(tmp.r)[1] <- sample.id.col
                                         tmp.r
                                     })
                        colnames(ret)[1] <- dataset.name.col
                        ret$id <- paste0(ret[, dataset.name.col], "-", ret[, sample.id.col])                        
                        ret
                    })
          })

## Calculate both pearson and spearman correlation over bootstraps, within dataset and cell type
bootstrapped.cors <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              flag <- res[, subchallenge.col] == sub.challenge
              tmp <- res[flag, ]
              flag <- !is.na(tmp[, measured.col])
              tmp <- tmp[flag, ]
              methods <- unique(tmp[, model.id.col])
              names(methods) <- methods
              indices <- 1:n.bootstraps
              names(indices) <- indices
              ret.i <-
                  ldply(indices,
                        .fun = function(i) {
                            ret.all <-
                                ldply(methods,
                                      .fun = function(method) {
                                          ret.method <-
                                              ldply(tbls.by.cell[[sub.challenge]],
                                                    .fun = function(df.in) {
                                                        flag <- df.in[, model.id.col] == method
                                                        df <- df.in[flag, ]
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
                                                                         mts <- c("pearson", "spearman")
                                                                         names(mts) <- mts
                                                                         vals <- llply(mts,
                                                                                       .fun = function(cor.method) {
                                                                                           if(any(is.na(df.ct[, prediction.col]))) {
                                                                                               if(!all(is.na(df.ct[, prediction.col]))) {
                                                                                                   stop("Some but not all NA\n")
                                                                                               }
                                                                                               return(NA)
                                                                                           }
                                                                                           if(var(as.numeric(df.ct[, prediction.col])) == 0) { return(0) }

                                                                                           val <- cor(as.numeric(df.ct[, prediction.col]),
                                                                                                      as.numeric(df.ct[, measured.col]),
                                                                                                      method = cor.method)
                                                                                           val
                                                                                       })
                                                                         as.data.frame(vals)
                                                                     })
                                                        colnames(tmp)[1] <- cell.type.col
                                                        tmp
                                                    })
                                          score
                                      })
                            colnames(ret.all)[1] <- model.id.col
                            ret.all
                        })
              colnames(ret.i)[1] <- "boot.i"
              ret.i
          })

bootstrapped.scores <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              df <- bootstrapped.cors[[sub.challenge]]
              un <- unique(res.latest[, c(model.id.col, method.name.col)])
              df <- merge(df, un)
              ## Average over cell type (within method and dataset and bootstrap sample)
              df <- ddply(df,
                          .variables = c(model.id.col, dataset.name.col, "boot.i"),
                          .fun = function(tmp) {
                              data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                          })
              ## Average over dataset (within method and bootstrap sample)
              df <- ddply(df,
                          .variables = c(model.id.col, "boot.i"),
                          .fun = function(tmp) {
                              data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                          })
              df
          })


mean.bootstrapped.scores <-
    llply(bootstrapped.scores,
          .fun = function(df) {
              ## Average over bootstraps
              df <- ddply(df,
                          .variables = c(model.id.col),
                          .fun = function(tmp) {
                              data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman))
                          })
              o <- order(df$pearson, decreasing = TRUE)
              df <- df[o,]
              df
          })

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))

for(sub.challenge in sub.challenges) {
    scores <- bootstrapped.scores[[sub.challenge]]
    mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
    un <- unique(res.latest[, c(model.id.col, method.name.col)])    
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

    png(paste0("validation-score-boxplots-", sub.challenge, ".png"), width = 2 * 480)
    title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
    grid.arrange(g1, g2, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))
    d <- dev.off()
}

calculate.empirical.bayes <-
    function(df, col.id, numerator.id, denominator.id, sample.id.cols, score.col) {
        flag <- df[, col.id] %in% c(numerator.id, denominator.id)
        tmp <- df[flag, ]
        n.tot <- nrow(tmp) / 2
        n.num <-
            sum(unlist(
                dlply(tmp, .variables = sample.id.cols,
                      .fun = function(df.2) {
                          if(nrow(df.2) != 2) { stop("Was expecting 2 rows\n") }
                          if(any(is.na(df.2[, score.col]))) { stop("Got NA scores\n") }
                          rownames(df.2) <- df.2[, col.id]
                          diff <- df.2[numerator.id, score.col] -
                              df.2[denominator.id, score.col]
                          if(diff > 0) { return(1) }
                          return(0)
                      })))
        n.num / (n.tot - n.num)
    }

bayes.factors <-
    llply(sub.challenges,
          .fun = function(sub.challenge) {
              mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
              mean.scores <- na.omit(mean.scores)
              scores <- bootstrapped.scores[[sub.challenge]]
              numerator.indx <- 1
              numerator.id <- mean.scores[numerator.indx, "modelId"]
              indices <- 1:nrow(mean.scores)
              indices <- indices[!(indices == numerator.indx)]
              names(indices) <- mean.scores[indices, "modelId"]
              ret <-
                  ldply(indices,
                        .fun = function(i) {
                            denominator.id <- mean.scores[i, "modelId"]
                            methods <- c("pearson", "spearman")
                            names(methods) <- methods
                            res <- llply(methods,
                                         .fun = function(method) {
                                             bf <- calculate.empirical.bayes(scores, col.id = "modelId",
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

stop("stop here\n")

source("../utils.R")

means.by.cell.type.method <-
    llply(bootstrapped.cors,
          .fun = function(df) {
              methods <- c("pearson", "spearman")
              na.rm <- FALSE
              names(methods) <- methods
              res <- llply(methods,
                           .fun = function(method) {
                               ## first, average over bootstrap
                               ret <- ddply(df, .variables = c(model.id.col, cell.type.col, dataset.name.col),
                                            .fun = function(df) {
                                                data.frame(cor = mean(df[, method], na.rm=na.rm))
                                            })
                               ## now, average over dataset
                               ret2 <- ddply(ret, .variables = c(model.id.col, cell.type.col),
                                             .fun = function(df) {
                                                 data.frame(cor = mean(df$cor, na.rm=na.rm))
                                             })
                           })
          })


for(sub.challenge in sub.challenges) {
    means <- means.by.cell.type.method[[sub.challenge]][["pearson"]]
    un <- unique(res.latest[, c(model.id.col, method.name.col)])    
    means <- merge(means, un)

    g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                            id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
    g <- g + ggtitle(paste0(firstup(sub.challenge), "-Grained Sub-Challenge"))

    png(paste0("validation-cell-heatmap-", sub.challenge, ".png"), width = 2 * 480)
    print(g)
    d <- dev.off()
}

save.image(".Rdata.bootstrap.validation")
