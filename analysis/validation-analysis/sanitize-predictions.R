suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

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

## Sanitize predictions from baseline (except cibersortX, which was
## run separately) and participant models. i.e., assign them sensible names
## and submission rounds.

## The predictions are those from the 'rerun-validation' / 'post-competitive' phase, i.e.,
## where the coarse- and fine-grained challenge use the same data.
## This is as opposed to the original competitive phase (against which the
## methods were ranked) where the coarse- and fine-grained challenges
## differed.

folder.synId <- "syn22320184"
predictions.synId <- "syn22314641"

## Read in the "post-competitive" predictions 
obj <- synGet(predictions.synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

## Read in the final predictions from the competitive phase (i.e., when coarse- and fine-grained
## datasets differed)
synId <- "syn22149603"
obj <- synGet(synId, downloadFile=TRUE)
res.comp <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

## Use the competitive phase results to translate objectIds to teams / submitters
## The final name in the "path" of the repo_name for the post-competitive phase is the
## objectId in the competitive phase (except for the baselines)
res.all$objectId <-
    unlist(llply(as.character(res.all$repo_name),
                 .fun = function(str) {
                     strs <- unlist(strsplit(str, "/"))
                     strs[length(strs)]
                 }))

res.all$comparator <-
    unlist(lapply(as.character(res.all$objectId),
                  function(str) {
                      comps <- get.comparators()
                      for(comp in comps) {
                          if(grepl(str, pattern=comp, ignore.case=TRUE)) { return(TRUE) }
                      }
                      return(FALSE)
                  }))

## Ensure that all objectIds match between competitive and post-competitive
flag <- grepl(res.comp$repo_name, pattern="baseline") | (res.comp$objectId %in% res.all$objectId)
if(!all(flag)) {
    stop("Some competitive objectIds are not in post-competitive results\n")
}

flag <- (res.all$comparator == TRUE) | (res.all$objectId %in% res.comp$objectId)
if(!all(flag)) {
    stop("Some post-competitive objectIds are not in competitive results\n")
}

safe.merge <- function(x, y, ...) {

    orig.nrow <- nrow(x)
    x <- merge(x, y, ...)
    new.nrow <- nrow(x)
    if(orig.nrow != new.nrow) {
        stop("Changed rows\n")
    }
    x
}

res.all <- safe.merge(res.all, unique(res.comp[, c("objectId", "submitterId")]), all.x=TRUE)

## Append CIBERSORTx results
## cibersortx.res.file <- "validation-csx-all-gene-predictions.tsv"
add.cibersortx.results <- TRUE
if(add.cibersortx.results) {
    cibersortx.res.synId <- "syn22329798"
    cibersortx.res.file <- synGet(cibersortx.res.synId, downloadFile=TRUE)$path
    cibersortx.res <- read.table(cibersortx.res.file, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
    cibersortx.res$objectId <- "CIBERSORTx"
    cibersortx.res$repo_name <- "CIBERSORTx"
    cibersortx.res$comparator <- TRUE
    cibersortx.res$submitterId <- NA
    measured.tbl <- unique(res.all[, c("dataset.name", "sample.id", "cell.type", "measured", "subchallenge")])
    cibersortx.res <- merge(cibersortx.res, measured.tbl)
    cibersortx.res <- cibersortx.res[, colnames(res.all)]
    res.all <- rbind(res.all, cibersortx.res)
}

## Assign the team name

team.name.tbl <- unique(res.all[, c("objectId", "subchallenge", "submitterId", "comparator")])
flag <- !(team.name.tbl$comparator == TRUE)

team.name.tbl$method.name <- NA
team.name.tbl[flag, "method.name"] <-
    unlist(lapply(as.character(team.name.tbl[flag, "submitterId"]),
                  function(str) translate.submitterId(str)))

team.name.tbl <-
    assign.baseline.names(team.name.tbl, from.col = "objectId", to.col = "method.name")


team.name.tbl <- simplify.submitter.names(team.name.tbl, col = "method.name")

## Assign the round
team.name.tbl <-
    assign.submission.rounds(team.name.tbl, object.id.col = "objectId",
                             context.cols = c("subchallenge", "submitterId"),
                             method.name.col = "method.name")

res.all <- merge(res.all, team.name.tbl, all.x = TRUE)

in.silico.validation.metadata.synId <- "syn22013519"
obj <- synGet(in.silico.validation.metadata.synId, downloadFile=TRUE)
in.silico.validation.metadata <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
in.silico.validation.metadata <- in.silico.validation.metadata[, c("tumor.type", "mixture.type", "dataset.name")]
colnames(in.silico.validation.metadata) <- c("tumor.type", "distribution.type", "dataset")
in.silico.validation.metadata$mixture.type <- "In Silico"

suppressPackageStartupMessages(p_load("openxlsx"))
in.vitro.validation.metadata <- get.validation.metadata()
in.vitro.validation.metadata <- unique(in.vitro.validation.metadata[, c(2,4,5)])
colnames(in.vitro.validation.metadata) <- c("tumor.type", "distribution.type", "dataset")
in.vitro.validation.metadata$mixture.type <- "In Vitro"

validation.metadata <- rbind(in.silico.validation.metadata, in.vitro.validation.metadata)

## Was this an in silico or an in vitro dataset?
flag <- res.all$dataset %in% validation.metadata$dataset
if(!all(flag)) {
    stop("Some datasets missing\n")
}

## Merge in distribution.type (random vs in biological), cancer type (tumor.type), and mixture.type (in silico vs in vitro)
res.all <- safe.merge(res.all, validation.metadata, all.x = TRUE, by.x = "dataset.name", by.y = "dataset")

## Store the results in synapse
file <- "rerun-validation-sanitized-predictions.csv"
write.table(file = file, res.all, col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
cat(paste0("Storing ", file, " to synapse\n"))
f <- File(file, parentId = folder.synId, synapseStore = TRUE)
synStore(f)

cat("Exiting successfully\n")
q(status=0)

stop("stop")

## TODO

## recreate: generate-in-silico-admixtures.R

## TODO

## Sensitivity and specificity
## - write combine sensitivity results
## - plot sensitivity heatmap based on combined sensitivity results (plot-limit-of-detection.R)

## Create plots
## - Fig 1: coarse-grained box and bar (rerun-validation-score-box-and-barplots-coarse-round-latest.png; done)
## - Fig 2a: coarse-grained heatmap (rerun-validation-bootstrap-cell-heatmap-coarse-round-latest.png; done)
## - Fig 2b: fine-grained heatmap (rerun-validation-bootstrap-cell-heatmap-fine-round-latest.png; done)
## - Fig 4: specificity heatmap (open all-genes-validation-spillover-all-coarse-grained.png; done)

## - Fig 3a: sensitivity example (validation-lod-ct-Random-B-cells-EPIC-coarse.png)
## - Fig 3b: sensitivity heatmap (validation-lod-summary-Random-coarse.png)



subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
model.id.col <- "modelId"
method.name.col <- "repo_name"
object.id.col <- "objectId"
submitter.id.col <- "submitterId"

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


cat(paste0("nrow before assign.rounds = ", nrow(res.all), "\n"))
context.cols <- c(subchallenge.col, submitter.id.col)
res.all <- assign.submission.rounds(res.all, object.id.col, context.cols, method.name.col)

print(table(res.all$submission))

my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

calculate.bootstraps <- function(res, n.bootstraps = 1000) {

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
                  llply(1:n.bootstraps, .parallel = TRUE,
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
    bootstraps
}

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

bootstraps <- calculate.bootstraps(res.all)

do.bootstrap.analysis <-
    function(res.input, boostraps, 
             method.name.col, object.id.col, submitter.id.col,
             model.id.col, subchallenge.col, measured.col, cell.type.col,
             dataset.name.col, sample.id.col, prediction.col,
             round.col, round = "latest",
             postfix) {

        submitter.tbl <- unique(res.input[, c(submitter.id.col, method.name.col, round.col, subchallenge.col), drop = FALSE])
        submitter.name.col <- "submitter"
        submitter.tbl[, submitter.name.col] <-
            unlist(lapply(submitter.tbl[, submitter.id.col],
                          function(x) translate.submitterId(x)))
        
        submitter.tbl <- simplify.submitter.names(submitter.tbl, col = submitter.name.col)
        submitter.tbl <- assign.baseline.names(submitter.tbl, from.col = method.name.col, to.col = submitter.name.col)

        lvls <- c("1", "2", "latest")
        submitter.tbl[, round.col] <- factor(as.character(submitter.tbl[, round.col]),
                                         levels = lvls)
        or <- order(submitter.tbl[, round.col])
        submitter.tbl <- submitter.tbl[or, ]
        submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
        flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
        submitter.tbl <- submitter.tbl[flag, ]
        flag <- !duplicated(submitter.tbl[, c(submitter.name.col, subchallenge.col)], fromLast = FALSE)
        submitter.tbl <- submitter.tbl[flag, ]
        
        res.round <- merge(res.input, submitter.tbl, by = c(submitter.id.col, method.name.col, subchallenge.col, round.col))
        
        ## baseline.method.flag <- define.baseline.method.flag(res.input, method.name.col)

        ## ## Always take the latest baseline method results
        ## non.baseline.method.flag <- !baseline.method.flag
        ## baseline.method.flag <- baseline.method.flag & (res.input[, round.col] == "latest")
        ## non.baseline.method.flag <- non.baseline.method.flag & (res.input[, round.col] == round)
        ## res.round <- res.input[baseline.method.flag | non.baseline.method.flag, ]
        
        res.round[, model.id.col] <- paste0(as.character(res.round[, method.name.col]), "-", as.character(res.round[, submitter.id.col]), "-", as.character(res.round[, object.id.col]))

        method.id.col <- model.id.col
        method.id.col <- submitter.name.col
        
        un <- unique(res.round[, c(method.name.col, submitter.id.col)])
        flag <- my.dup(un[, method.name.col])
        dup.repo.names <- un[flag, method.name.col]
        flag <- res.round[, method.name.col] %in% dup.repo.names
        res.round[flag, method.name.col] <- paste0(as.character(res.round[flag, method.name.col]), "-", as.character(res.round[flag, object.id.col]))
        
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
                      names(methods) <- methods
                      n.bootstraps <- length(bootstraps[[sub.challenge]])
                      indices <- 1:n.bootstraps
                      names(indices) <- indices
                      ret.i <-
                          ldply(indices, .parallel = TRUE,
                                .fun = function(i) {
                                    ret.all <-
                                        ldply(methods,
                                              .fun = function(method) {
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
                                    colnames(ret.all)[1] <- method.id.col
                                    ret.all
                                })
                      colnames(ret.i)[1] <- "boot.i"
                      ret.i
                  })

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
        

        cat(paste0("Calculating boxplots\n"))
        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            mean.scores <- mean.bootstrapped.scores[[sub.challenge]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, unique(c(method.id.col, method.name.col, submitter.name.col))])
            mean.scores <- merge(mean.scores, un)
            scores <- merge(scores, un)    
            o <- order(mean.scores$pearson)
            mean.scores <- mean.scores[o, ]
            ## scores[, method.name.col] <- factor(scores[, method.name.col], levels = mean.scores[, method.name.col])
            scores[, submitter.name.col] <- factor(scores[, submitter.name.col], levels = mean.scores[, submitter.name.col])
            scores <- na.omit(scores)
            
            g1 <- ggplot(data = scores)
            g1 <- g1 + geom_boxplot(aes_string(x = submitter.name.col, y = "pearson"))
            g1 <- g1 + coord_flip()
            g1 <- g1 + xlab("Method")
            g1 <- g1 + ylab("Pearson Correlation")
            g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))
            
            g2 <- ggplot(data = scores)
            g2 <- g2 + geom_boxplot(aes_string(x = submitter.name.col, y = "spearman"))
            g2 <- g2 + coord_flip()
            g2 <- g2 + xlab("Method")
            g2 <- g2 + ylab("Spearman Correlation")
            g2 <- g2 + theme(text = element_text(size=18))    
            
            png(paste0("validation-score-boxplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
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

        cat(paste0("Calculating mean by cell type\n"))        
        means.by.cell.type.method <-
            llply(bootstrapped.cors,
                  .fun = function(df) {
                      methods <- c("pearson", "spearman")
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
        
        cat(paste0("Plotting heatmaps\n"))
        for(sub.challenge in sub.challenges) {
            means <- means.by.cell.type.method[[sub.challenge]][["pearson"]]
            flag <- res.round[, subchallenge.col] == sub.challenge            
            un <- unique(res.round[flag, unique(c(method.id.col, method.name.col, submitter.name.col))])    
            means <- merge(means, un)
            
            g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                    id.var = submitter.name.col, cell.type.var = cell.type.col, cor.var = "cor")
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
            png(paste0("validation-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
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
for(round in c("1", "2", "3", "latest")) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))
    results[[round]] <- do.bootstrap.analysis(res.all, bootstraps, 
                                              method.name.col, object.id.col, submitter.id.col,
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
        file <- paste0("bayes-factors-round-", round, "-sc-", sub.challenge, "-vs-",
                       make.names(best.team), ".tsv")
        write.table(file = file, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    }
}

save.image(".Rdata.bootstrap.validation")

cat("Exiting successfully\n")
q(status=0)
