suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))
suppressPackageStartupMessages(p_load("reshape2"))

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

if(FALSE) {
        # Add ensemble here -- no need. We report bootstrapped scores
        # > colnames(results[["3"]][["res.round"]])
        # [1] "method.name"       "subchallenge"      "submission"       
        # [4] "dataset.name"      "sample.id"         "cell.type"        
        # [7] "objectId"          "comparator"        "submitterId"      
        # [10] "repo_name"         "tumor.type"        "distribution.type"
        # [13] "mixture.type"      "measured"          "prediction"       
        res.ensemble <-
          ddply(res.round[, !(colnames(res.round) %in% c("objectId", "comparator", "submittedId", "repo_name"))],
                .variables = c(round.col, subchallenge.col, dataset.name.col),
                .fun = function(df.in) {
                         tres <- 
                           ddply(df.in, .variables = method.name.col,
                                 .fun = function(df) {
                                          df.ret <- data.frame(prediction = df[, prediction.col],
                                                               pred.rank = rank(df[, prediction.col], ties.method="first"),
                                                               sample.id = df[, sample.id.col])
                                          colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                                          df.ret
                                        })
                         # take the consensus (here, just mean) rank across methods
                         cons <- 
                           ddply(tres, .variables = sample.id.col,
                                 .fun = function(df) {
                                          df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                                          colnames(df.ret)[1] <- prediction.col
                                          df.ret
                                        })
                         tmp <- unique(df.in[, c(sample.id.col, measured.col, subchallenge.col, round.col, dataset.name.col, cell.type.col, "tumor.type", "distribution.type", "mixture.type")])
                         cons <- merge(cons, tmp)
                         cons$method.name <- "ensemble"
                         cons$objectId <- "dummy"
                         cons$comparator <- FALSE
                         cons$submitterId <- "bwhite"
                         cons$repo_name <- "dummy"
                         ret <- cons[, c(method.name.col, subchallenge.col, round.col, dataset.name.col, sample.id.col, cell.type.col, 
                                         "objectId", "comparator", "submitterId", "repo_name", "tumor.type", "distribution.type", "mixture.type",
                                         measured.col, prediction.col)]
                         ret
                       })
        print(head(res.round))
        print(head(res.ensmble))
        res.round <- rbind(res.round, res.ensemble)
} #end if(FALSE)
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
        cat(paste0("Calculating bootstrap correlations (ensemble)\n"))
        bootstrapped.cors.ensemble <-
            llply(sub.challenges, 
                  .fun = function(sub.challenge) {
                      n.bootstraps <- length(bootstraps[[sub.challenge]])
                      indices <- 1:n.bootstraps
# indices <- 1:10
                      names(indices) <- indices

                      ret.i <-
                        ldply(indices,
                              .fun = function(i) {
                                       ret.cons <-
                                         ddply(tbls[[sub.challenge]], .variables = c(cell.type.col, dataset.name.col),
                                               .fun = function(df.in) {
                                                        # method.name is row
                                                        x <- acast(df.in[, c(prediction.col, method.name.col, sample.id.col)], 
                                                                   formula = paste0(method.name.col, " ~ ", sample.id.col), value.var = prediction.col)

                                                        colnames(x) <- paste0(df.in[1, dataset.name.col], "-", colnames(x))
                                                        sample.ids <- bootstraps[[sub.challenge]][[i]]$id
                                                        sample.ids <- sample.ids[sample.ids %in% colnames(x)]
                                                        x <- x[, sample.ids]
                                                        m <- melt(x)
                                                        colnames(m) <- c(method.name.col, sample.id.col, prediction.col)
                                                        tres <- 
                                                          ddply(m, .variables = method.name.col,
                                                                .fun = function(df) {
                                                                         df.ret <- data.frame(prediction = df[, prediction.col],
                                                                                              pred.rank = rank(df[, prediction.col], ties.method="first"),
                                                                                              sample.id = df[, sample.id.col])
                                                                         colnames(df.ret) <- c(prediction.col, "pred.rank", sample.id.col)
                                                                         df.ret
                                                                       })
                                                        # take the consensus (here, just mean) rank across methods
                                                        cons <- 
                                                          ddply(tres, .variables = sample.id.col,
                                                                .fun = function(df) {
                                                                         df.ret <- data.frame(cons.rank = mean(df$pred.rank))
                                                                         df.ret
                                                                       })
                                                        tmp <- unique(df.in[, c(sample.id.col, measured.col)])
                                                        tmp[, sample.id.col] <- paste0(df.in[1, dataset.name.col], "-", tmp[, sample.id.col])
                                                        cons <- merge(cons, tmp)
                                                        ret <- data.frame(method.name = "ensemble", pearson = NA, rmse = NA, spearman = cor(cons$cons.rank, cons[, measured.col]))
                                                        colnames(ret) <- c(method.name.col, "pearson", "rmse", "spearman")
                                                        ret
                                                      })
                                     })
        # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" 
                                   colnames(ret.i)[1] <- "boot.i"
                                   ret.i <- ret.i[, c(method.name.col, "boot.i", dataset.name.col, cell.type.col, "pearson", "spearman", "rmse")]
                                   ret.i 
                            })
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

        for(nm in names(bootstrapped.cors)) {
          bootstrapped.cors[[nm]] <- rbind(bootstrapped.cors[[nm]], bootstrapped.cors.ensemble[[nm]])
        }

        print(colnames(bootstrapped.cors[[1]]))
        # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" 
        print(head(bootstrapped.cors[[1]]))

        cat(paste0("Calculating bootstrapped scores\n"))
        bootstrapped.scores <-
            llply(sub.challenges, .parallel = TRUE,
                  .fun = function(sub.challenge) {
                      df <- bootstrapped.cors[[sub.challenge]]
                      flag <- res.round[, subchallenge.col] == sub.challenge
                      un <- unique(res.round[flag, unique(c(method.id.col, method.name.col))])
                      df <- merge(df, un)
                      ## Average over cell type (within method and dataset and bootstrap sample)
                      df <- ddply(df,
                                  .variables = c(method.id.col, dataset.name.col, "boot.i"),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse))
                                  })
                      ## Average over dataset (within method and bootstrap sample)
                      df <- ddply(df,
                                  .variables = c(method.id.col, "boot.i"),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse))
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
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse))
                                  })
                      o <- order(df$pearson, decreasing = TRUE)
                      df <- df[o,]
                      df
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

        means.over.bootstrap <-
            llply(bootstrapped.cors,
                  .fun = function(df) {
                      methods <- c("pearson", "spearman", "rmse")
                      na.rm <- FALSE
                      names(methods) <- methods
                      res <- llply(methods,
                                   .fun = function(method) {
                                       ## average over dataset
                                       ret <- ddply(df, .variables = c(method.id.col, cell.type.col, dataset.name.col),
                                                    .fun = function(df) {
                                                        data.frame(cor = mean(df[, method], na.rm=na.rm))
                                                    })
                                   })
                  })
        
        cat(paste0("Calculating mean by cell type\n"))
        if(FALSE) {
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
        }

        means.by.cell.type.method <-
                llply(means.over.dataset,
                      .fun = function(df) {
                          methods <- c("pearson", "spearman", "rmse")
                          na.rm <- FALSE
                          names(methods) <- methods
                          res <- llply(methods,
                                       .fun = function(method) {
                                           ## average over bootstrap (means.over.dataset has already been averaged over dataset)
                                           ret <- ddply(df[[method]], .variables = c(method.id.col, cell.type.col),
                                                        .fun = function(df) {
                                                            data.frame(cor = mean(df$cor, na.rm=na.rm))
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

        ret.list <- list("res.round" = res.round,
                         "bootstrapped.cors" = bootstrapped.cors,
                         "bootstrapped.scores" = bootstrapped.scores,
                         "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "means.by.cell.type.method" = means.by.cell.type.method,
                         "means.over.dataset" = means.over.dataset,
                         "means.over.bootstrap" = means.over.bootstrap,                         
                         "top.performers" = top.performers,
                         "bayes.factors" = bayes.factors)
                         
        return(ret.list)
        
}

results <- list()
## for(round in c("1", "2", "3", "latest")) {
for(round in c("1", "2", "3", "latest")) {
# for(round in c("1")) {
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

## save.image(".Rdata.plot.bootstrap")

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

cat("Exiting successfully\n")
q(status=0)


