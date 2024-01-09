load.challenge.validation.data <- function() { 
  # Location of Challenge validation admixtures on Synapse
  validation.data.synpase.folder.id <- "syn21821096"

  # List all files in the Challenge validation folder
  children <- synGetChildren(validation.data.synpase.folder.id)
  l <- as.list(children)
  df <- do.call(rbind.data.frame, l)

  # Load the "input.csv" file that lists the expression matrices in each
  # of the Challenge validation datasets
  input.file.synId <- as.character(subset(df, name == "input.csv")[, "id"])
  obj <- synGet(input.file.synId, downloadFile=TRUE)
  input.tbl <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

  # Subset to the expression data with Hugo symbols
  input.files <- input.tbl[, c("dataset.name", "hugo.expr.file")]
  df <- merge(df, input.files, by.x = "name", by.y = "hugo.expr.file")

  # Read in each of the datasets
  synIds <- as.character(df$id)
  names(synIds) <- df$dataset.name

  # Download each dataset from Synapse and read it in
  datasets <- 
    llply(synIds, 
          .fun = function(synId) {
                   obj <- synGet(synId, downloadFile=TRUE)
                   mat <- as.data.frame(fread(obj$path))
                   rownames(mat) <- mat$Gene
                   mat <- subset(mat, select = -Gene)
                 })

  datasets
}

# this defines methods.to.exclude
source("methods-to-exclude.R")

# Return a list with named entries:
# predictions: the predictions from the comparator and participant methods, with columns
#              dataset.name, subchallenge, sample.id, cell.type, prediction, method.name, and submission (1, 2, 3, or latest)
# dataset.anno: the description of each dataset, with columns
#              dataset.name, tumor.type, distribution.type, mixture.type
# ground.truth: the ground truth, with columns:
#              dataset.name, subchallenge, sample.id, cell.type, measured
# In all cases, limit results to target.submissions and target.subchallenges
load.challenge.validation.results <- function(target.submissions = c("1", "2", "3", "latest"), target.subchallenges = c("fine", "coarse"),
                                              method.name.col = "method.name") {
  synId <- "syn22320329"
  obj <- synGet(synId, downloadFile=TRUE)
  res.all <- as.data.frame(fread((obj$path)))

  flag <- res.all[,method.name.col] %in% methods.to.exclude
  cat(paste0("Excluding methods: ", paste(unique(res.all[flag, method.name.col]), collapse = ", "), "\n"))
  res.all <- res.all[!flag,]

  res.all <- subset(res.all, submission %in% target.submissions)
  res.all <- subset(res.all, subchallenge %in% target.subchallenges)

  preds <- res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "prediction", "method.name", "submission")]
  anno <- unique(res.all[, c("dataset.name", "tumor.type", "distribution.type", "mixture.type")])
  rownames(anno) <- NULL
  gt <- unique(res.all[, c("dataset.name", "subchallenge", "sample.id", "cell.type", "measured")])
  rownames(gt) <- NULL

  ret <- list("predictions" = preds, "dataset.anno" = anno, "ground.truth" = gt)
  ret
}

## Ensure we have a prediction (even if it is NA) for all cell types in all datasets by all methods
make.missing.predictions.na <- function(res.all, cell.type.col = "cell.type", prediction.col = "prediction",
                                          measured.col = "measured", subchallenge.col = "subchallenge",
                                          sample.id.col = "sample.id", dataset.name.col = "dataset.name",
                                          round.col = "submission") {
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
  res.all
}

# Aggregate predictions from raw cell types as output by a deconvolution method
# into those cell types expected of the Challenge.
# Aggregation is performed via summing.
# The translation between raw and Challenge subtypes is encoded in the translation table,
# where the raw.cell.type column gives the raw cell type name and cell.type gives the Challenge cell type.
# Results in res_df are described by columns method.name, dataset.name, sample.id, cell.type, and prediction
aggregate.cell.type.predictions <- function(res_df, translation_df) {
  nm <- unique(res_df$method.name)
  res_df <- res_df %>%
    dplyr::inner_join(translation_df) %>%
    dplyr::select(dataset.name, sample.id, cell.type, prediction) %>%
    dplyr::group_by(dataset.name, sample.id, cell.type) %>%
    dplyr::summarise(prediction = sum(prediction)) %>%
    dplyr::ungroup()
  res_df <- as.data.frame(res_df)
  res_df <- cbind(method.name = nm, res_df)
}

sub.challenges <- list("coarse" = "coarse", "fine" = "fine")

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
             subchallenge.col, measured.col, cell.type.col,
             dataset.name.col, sample.id.col, prediction.col,
             round.col, round = "latest") {

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
                  .parallel = FALSE,
                  .fun = function(sub.challenge) {
                      n.bootstraps <- length(bootstraps[[sub.challenge]])
                      indices <- 1:n.bootstraps
# indices <- 1:10
                      names(indices) <- indices

                      ret.i <-
                        ldply(indices,
                              .parallel = TRUE,
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
                                                        ret <- data.frame(method.name = "ensemble", pearson = NA, rmse = NA, spearman = cor(cons$cons.rank, cons[, measured.col]), pearson.fc = NA)
                                                        colnames(ret) <- c(method.name.col, "pearson", "rmse", "spearman", "pearson.fc")
                                                        ret
                                                      })
                                     })
        # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
                                   colnames(ret.i)[1] <- "boot.i"
                                   ret.i <- ret.i[, c(method.name.col, "boot.i", dataset.name.col, cell.type.col, "pearson", "spearman", "rmse", "pearson.fc")]
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
                                                                                 mts <- c("pearson", "spearman", "rmse", "pearson.fc")
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
                                                                                                   } else if(comp.method == "pearson.fc") {
                                                                                                       # calculate the pearson correlation of log fold changes (across successive samples)
                                                                                                       # 1. Order the samples based on the their ground truth values g1 < g2 < … < gn (i.e., measured)
                                                                                                       # 2. Compute ground truth fold changes g2/g1, g3/g2, …, gn/gn-1
                                                                                                       # 3. Compute prediction fold changes p2/p1, p3/p2, …, pn/pn-1
                                                                                                       # 4. Compute the pearson correlation of these two vectors
                                                                                                       o <- order(measured, decreasing=FALSE)
                                                                                                       eps <- 10^-5
                                                                                                       gt.ordered <- measured[o] + eps
                                                                                                       # NB: the ordering is established by the _ground truth_ values and applied to those
                                                                                                       # values and the predicted values
                                                                                                       pred.ordered <- pred[o] + eps
                                                                                                       gt.fc <- unlist(lapply(2:length(gt.ordered), function(indx) gt.ordered[indx]/gt.ordered[indx-1]))
                                                                                                       pred.fc <- unlist(lapply(2:length(pred.ordered), function(indx) pred.ordered[indx]/pred.ordered[indx-1]))
                                                                                                       val <- cor(pred.fc, gt.fc, method = "pearson")
                                                                                                   } else if(comp.method == "rmse") {
                                                                                                       val <- sqrt(mean((pred-measured)^2))
                                                                                                   } else {
                                                                                                       stop(paste0("Unknown method ", comp.method, "\n"))
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
        # "method.name"  "boot.i"       "dataset.name" "cell.type"    "pearson"   "spearman"     "rmse" "pearson.fc"
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
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                                  })
                      ## Average over dataset (within method and bootstrap sample)
                      df <- ddply(df,
                                  .variables = c(method.id.col, "boot.i"),
                                  .fun = function(tmp) {
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
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
                                      data.frame(pearson = mean(tmp$pearson), spearman = mean(tmp$spearman), rmse = mean(tmp$rmse), pearson.fc = mean(tmp$pearson.fc))
                                  })
                      o <- order(df$pearson, decreasing = TRUE)
                      df <- df[o,]
                      df
                  })

        means.over.dataset <-
            llply(bootstrapped.cors,
                  .fun = function(df) {
                      methods <- c("pearson", "spearman", "rmse", "pearson.fc")
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
                      methods <- c("pearson", "spearman", "rmse", "pearson.fc")
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
                          methods <- c("pearson", "spearman", "rmse", "pearson.fc")
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
                          methods <- c("pearson", "spearman", "rmse", "pearson.fc")
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
    function(res, bootstrapped.scores, mean.boostrapped.scores, median.bootstrapped.scores,
             means.by.cell.type.method,
             means.over.dataset, method.anno.round,
             postfix, plot.spearman.distribution = FALSE) {

        top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
        comparator.methods <- get.comparators.cap()
	print(comparator.methods) 
        priority.methods <- unique(c(top.performers, comparator.methods))

        for(sub.challenge in sub.challenges) {
            print(sub.challenge)
            print(head(method.anno.round))
	    print(method.anno.round)
            print(subchallenge.col)
            flag <- is.na(method.anno.round[, subchallenge.col]) | ( method.anno.round[, subchallenge.col] == sub.challenge )
            print(flag)
            method.anno.round.sc <- method.anno.round[flag, c(method.name.col, "Output", "Method")]
print(method.anno.round.sc)
            bootstrapped.scores[[sub.challenge]] <-
                merge(bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
            median.bootstrapped.scores[[sub.challenge]] <-
                merge(median.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
            mean.bootstrapped.scores[[sub.challenge]] <-
                merge(mean.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
        }
        
        barplots <- list()
        for(sub.challenge in sub.challenges) {
            scores <- bootstrapped.scores[[sub.challenge]]
            median.scores <- median.bootstrapped.scores[[sub.challenge]]

            median.scores[, method.name.col] <- as.character(median.scores[, method.name.col])
            flag <- median.scores[, method.name.col] == "ensemble"
            median.scores[flag, method.name.col] <- "consensus rank"

            scores[, method.name.col] <- as.character(scores[, method.name.col])
            flag <- scores[, method.name.col] == "ensemble"
            scores[flag, method.name.col] <- "consensus rank"


            o <- order(median.scores$pearson)
            median.scores <- median.scores[o, ]
            flag <- ( !is.na(scores$pearson) & !is.na(scores$spearman) ) | ( scores[,method.name.col] == "consensus rank")
            scores <- scores[flag, ] 
            scores[, method.name.col] <- factor(scores[, method.name.col], levels = median.scores[, method.name.col])

            lvls <- levels(scores[,method.name.col]) 
            lvls <- lvls[lvls %in% scores[, method.name.col]]
            bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
cat("scores[, method.name.col]\n")
print(lvls)
print(bold.labels)
print(table(bold.labels))
#            scores <- na.omit(scores)
            flag <- ( !is.na(median.scores$pearson) & !is.na(median.scores$spearman) ) | ( median.scores[,method.name.col] == "consensus rank") | ( median.scores[,method.name.col] == "ensemble")
            median.scores[, method.name.col] <- factor(median.scores[, method.name.col], levels = median.scores[, method.name.col])
#            median.scores <- na.omit(median.scores)
            median.scores <- median.scores[flag, ] 
            

            g1 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson"))
            g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
            g1 <- g1 + coord_flip()
            g1 <- g1 + xlab("Method")
            ## g1 <- g1 + ylab("Pearson Correlation")
            g1 <- g1 + ylab("Pearson")
            g1 <- g1 + ylim(c(-0.25, 1))
            g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
            # g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))
            # g1 <- g1 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
            # g1 <- g1 + theme(axis.text.y=element_markdown())

            if(plot.spearman.distribution) {
              g2 <- ggplot(data = scores, aes_string(x = method.name.col, y = "spearman"))
              g2 <- g2 + geom_boxplotMod(fill = "#E69F00")
            } else{
              g2 <- ggplot(data = median.scores)
              g2 <- g2 + geom_col(aes_string(x = method.name.col, y = "spearman"), fill = "#E69F00")
            }
            g2 <- g2 + coord_flip()
            g2 <- g2 + xlab("Method")
            ## g2 <- g2 + ylab("Spearman Correlation")
            g2 <- g2 + ylab("Spearman")
            g2 <- g2 + ylim(c(-0.25, 1))            
            g2 <- g2 + theme(text = element_text(size=18))    
            g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.ticks.y = element_blank())

            g3 <- ggplot(data = scores, aes_string(x = method.name.col, y = "pearson.fc"))
            g3 <- g3 + geom_boxplotMod(fill = "#56B4E9")
            g3 <- g3 + coord_flip()
            g3 <- g3 + xlab("Method")
            ## g3 <- g3 + ylab("Pearson Correlation")
            g3 <- g3 + ylab("Pearson (Fold Change)")
            g3 <- g3 + ylim(c(-0.25, 1))
            g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20), axis.text.y = element_text(face = bold.labels))
            # g3 <- g3 + theme(text = element_text(size=18), title = element_text(size = 20))
            # g3 <- g3 + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
            # g3 <- g3 + theme(axis.text.y=element_markdown())
            
     
            tmp <- scores[, c(method.name.col, "Output", "Method")]            
            ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))

            full.plot <- ret[["full.plot"]]
            for.first.legend <- ret[["legends"]][["Method"]]
            for.second.legend <- ret[["legends"]][["Output"]]
            
            legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
            
            ## pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,1,0.5,0.5))
            
            barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
            barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            
            barplots[[paste0(sub.challenge,  "-pearson.fc")]] <- g3            
            barplots[[paste0(sub.challenge,  "-anno")]] <- full.plot
            barplots[[paste0(sub.challenge,  "-legend")]] <- legs

            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)            
            title <- paste0(title, " (", round.text, ")")
            ## png(paste0(figs.dir, "rerun-validation-score-box-and-barplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## grid.arrange(g1, g2, nrow=1, widths = c(3, 1), top = textGrob(title, gp = gpar(fontsize = 25)))
            ## d <- dev.off()
        }

##        ret.list <- list(
##            "barplots" = barplots)
##        return(ret.list)
                    
        
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

        cat(paste0("Plotting heatmaps\n"))
        heatmaps <- list()
	spearman.heatmaps <- list()
	pearson.fc.heatmaps <- list()
        method.levels <- list()
        cell.type.levels <- list()
        for(sub.challenge in sub.challenges) {
	for(cor.type in c("pearson", "spearman", "pearson.fc")) {
            means <- means.by.cell.type.method[[sub.challenge]][[cor.type]]

            method.levels[[paste0(sub.challenge, "-", cor.type)]] <-
                calculate.method.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
            cell.type.levels[[paste0(sub.challenge, "-", cor.type)]] <-
                calculate.cell.type.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
            
            exclude.method <-
              ddply(means, .variables = c(method.name.col),
                    .fun = function(df) {
                             exclude <- all(is.na(df$cor))
                             data.frame(exclude = exclude)
                           })
            if(any(exclude.method$exclude)) {
              means <- means[!(means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
              method.levels[[paste0(sub.challenge, "-", cor.type)]] <- 
                method.levels[[paste0(sub.challenge, "-", cor.type)]][!(method.levels[[paste0(sub.challenge, "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
            }
            cat("means\n")
print(unique(means[, method.name.col]))
cat("levels\n")
print(method.levels[[paste0(sub.challenge, "-", cor.type)]])
            cor.type.label <- "Pearson\nCorrelation"
	    if(cor.type == "spearman") {
	      cor.type.label <- "Spearman\nCorrelation"
	    }
	    if(cor.type == "pearson.fc") {
	      cor.type.label <- "Pearson Correlation\n(Fold Change)"
	    }
            g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                    id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                                    second.col.summary.fun = "mean",
						    cor.type.label = cor.type.label,
                                                    method.levels = method.levels[[paste0(sub.challenge, "-", cor.type)]],
                                                    cell.type.levels = cell.type.levels[[paste0(sub.challenge, "-", cor.type)]], ids.to.bold = comparator.methods)
##            g <- plot.cell.type.correlation.strip.plots(means, show.corr.text = TRUE, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
            if(cor.type == "pearson") { 
              heatmaps[[sub.challenge]] <- g
            } else if(cor.type == "spearman") {
              spearman.heatmaps[[sub.challenge]] <- g	    
	    } else if(cor.type == "pearson.fc") {
              pearson.fc.heatmaps[[sub.challenge]] <- g
            } else {
              stop(paste0("Unknown cor.type ", cor.type, "\n"))
            }
            title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
            round.text <- make.round.text(round)
            title <- paste0(title, " (", round.text, ")")
            
            g <- g + ggtitle(title)
            ## png(paste0(figs.dir, "rerun-validation-bootstrap-cell-heatmap-", sub.challenge, postfix, ".png"), width = 2 * 480)
            ## print(g)
            ## d <- dev.off()
        } # for cor.type
	}

	for(cor.type in c("pearson", "spearman", "pearson.fc")) {
        coarse.means <- means.by.cell.type.method[["coarse"]][[cor.type]]
        fine.means <- means.by.cell.type.method[["fine"]][[cor.type]]        
        all.means <- rbind(coarse.means, fine.means)

        all.means <-
            ddply(all.means, .variables = c(method.name.col, cell.type.col),
                  .fun = function(df) {
                      data.frame(cor = summary.fun(df$cor))
                  })

        method.levels[[paste0("merged", "-", cor.type)]] <-
            calculate.method.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")

        cell.type.levels[[paste0("merged", "-", cor.type)]] <-
            calculate.cell.type.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        
        exclude.method <-
          ddply(all.means, .variables = c(method.name.col),
                .fun = function(df) {
                         exclude <- all(is.na(df$cor))
                         data.frame(exclude = exclude)
                       })
cat(paste0("Exclude: ", cor.type, "\n"))
print(exclude.method)
        if(any(exclude.method$exclude)) {
          all.means <- all.means[!(all.means[, method.name.col] %in% subset(exclude.method, exclude == TRUE)[, method.name.col]),]
          method.levels[[paste0("merged", "-", cor.type)]] <- 
            method.levels[[paste0("merged", "-", cor.type)]][!(method.levels[[paste0("merged", "-", cor.type)]] %in% subset(exclude.method, exclude == TRUE)[, method.name.col])]
        }

        cor.type.label <- "Pearson\nCorrelation"
	if(cor.type == "spearman") {
	  cor.type.label <- "Spearman\nCorrelation"
	}
	if(cor.type == "pearson.fc") {
	  cor.type.label <- "Pearson Correlation\n(Fold Change)"
	}

        g <- plot.cell.type.correlation.heatmap(all.means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
						cor.type.label = cor.type.label,						
                                                second.col.summary.fun = "mean",
                                                method.levels = method.levels[[paste0("merged", "-", cor.type)]],
                                                cell.type.levels = cell.type.levels[[paste0("merged", "-", cor.type)]], ids.to.bold = comparator.methods)
        merged.all.means <- all.means
	if(cor.type == "pearson") { 
          heatmaps[["merged"]] <- g
        } else if(cor.type == "spearman") {
          spearman.heatmaps[["merged"]] <- g	
	} else if(cor.type == "pearson.fc") {
          pearson.fc.heatmaps[["merged"]] <- g
        } else {
          stop(paste0("Unknown cor.type ", cor.type, "\n"))
        }
	} # for cor.type
	
        cat("Creating merged means.over.dataset\n")
        coarse.means <- means.over.dataset[["coarse"]][["pearson"]]
        fine.means <- means.over.dataset[["fine"]][["pearson"]]        
        all.means <- rbind(coarse.means, fine.means)

        all.means <-
            ddply(all.means, .variables = c(method.name.col, cell.type.col, "boot.i"),
                  .fun = function(df) {
                      data.frame(cor = summary.fun(df$cor))
                  })

        cat(paste0("Plotting strip plots\n"))
        nms <- list("coarse" = "coarse", "fine" = "fine", "coarse-priority" = "coarse-priority", "fine-priority" = "fine-priority",
                    "merged" = "merged", "merged-priority" = "merged-priority")
        strip.plots <-
            llply(nms,
                  .parallel = FALSE,
                  .fun = function(nm) {
                      sub.challenge <- NA
                      df <- NULL
                      entry <- NULL
                      lvl.entry <- NULL
                      if(grepl(nm, pattern="coarse")) {
                          entry <- "coarse"
                          lvl.entry <- paste0(entry, "-pearson")
                          df <- means.over.dataset[[entry]][["pearson"]]
                      }
                      if(grepl(nm, pattern="fine")) {
                          entry <- "fine"
                          lvl.entry <- paste0(entry, "-pearson")
                          df <- means.over.dataset[[entry]][["pearson"]]                          
                      }
                      if(grepl(nm, pattern="merged")) {
                          entry <- "merged"
                          lvl.entry <- paste0(entry, "-pearson")
                          df <- all.means
                      }

                      g <- NULL
                      if(grepl(nm, pattern="priority")) {
                          flag <- df[, method.name.col] %in% priority.methods
                          ret <- plot.strip.plots(df[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                                method.levels = method.levels[[lvl.entry]],
                                                cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                                label = "Pearson Correlation")
                          g <- ret[["g"]]
                          df <- ret[["df"]]
                          lvls <- levels(df[,method.name.col])
                          lvls <- lvls[lvls %in% df[,method.name.col]]
                          y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
print("here\n")
print(comparator.methods)
print(priority.methods)
print(lvls)
print(y.bold.labels)
                      } else {
                          # Exclude ensemble, which has NAs for pearson
                          flag <- !(df[, method.name.col] %in% c("consensus rank", "ensemble"))
                          ret <- plot.strip.plots(df[flag,], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                                method.levels = method.levels[[lvl.entry]],
                                                cell.type.levels = rev(cell.type.levels[[lvl.entry]]),
                                                label = "Pearson Correlation")
                          g <- ret[["g"]]
                          df <- ret[["df"]]
                          lvls <- levels(df[,method.name.col])
                          lvls <- lvls[lvls %in% df[,method.name.col]]
                          y.bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
print("there\n")
print(comparator.methods)
print(lvls)
print(y.bold.labels)


                      }
		      g <- g + theme(axis.text.y = element_text(face = y.bold.labels))
		      # g <- g + scale_y_discrete(labels = function(x) bold.highlight(x, comparator.methods))
                      # g <- g + theme(axis.text.y=element_markdown())
                  })

        ## "boxplots" = boxplots,
        ret.list <- list("median.bootstrapped.scores" = median.bootstrapped.scores,
                         "mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "barplots" = barplots,
                         "strip.plots" = strip.plots,
                         "heatmaps" = heatmaps,
                         "spearman.heatmaps" = spearman.heatmaps,			 
                         "pearson.fc.heatmaps" = pearson.fc.heatmaps,			 
			 "merged.all.means" = merged.all.means)			 

        ret.list
        
}


