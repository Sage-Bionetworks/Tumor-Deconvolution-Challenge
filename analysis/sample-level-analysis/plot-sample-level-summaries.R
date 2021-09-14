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

suppressPackageStartupMessages(p_load("xlsx"))

source("../utils.R")

set.seed(1234)

synLogin()

figs.dir <- "figs/"
dir.create(figs.dir, showWarnings = FALSE)

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

## Read in results from updated CIBERSORTx code that collapses sub-popuations (e.g., mem CD8 T and naive CD8 T)
## into their parental population (e.g., CD8 T) by correctly using sum rather than mean (as was used
## during the challenge). We will apply these changes to the CIBERSORT results (without re-running them).
## These files have columns cell.type, subchallenge, method.name, and revised.orig.ratio.
## method.name = CIBERSORTx (not CIBERSORT!) since CIBERSORTx combines the same sub-populations into parental units
## and was actually re-run. 
## We should revise the old cs result by multiplying it be revised.orig.ratio (i.e., convert mean to sum)
cs.revised.synIds <- list("coarse" = "syn26141656", "fine" = "syn26141634")
## syn26141656: specificity-coarse-csx-all-gene-predictions-revised-orig-ratios.tsv
## syn26141634: specificity-fine-csx-all-gene-predictions-revised-orig-ratios.tsv
## Hmmm ... for some reason, both of these files contain both coarse and fine-grained mappings.
## Ahh ... the different files result from the coarse- vs fine-grained subchallenges / datasets.
## But, CIBERSORTx was run so as to make coarse- or fine-grained predictions for both challenges.
## Proably the files are identical, but to be sure, only take the fine-grained mappings from the
## fine-grained challenge, etc.
revised.dfs <- 
  ldply(cs.revised.synIds, 
        .fun = function(synId) { 
                 obj <- synGet(synId, downloadFile=TRUE)
                 df <- read.table(obj$path, sep="\t", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
               })
colnames(revised.dfs)[1] <- "subchallenge.run"
revised.dfs$method.name <- "CIBERSORT"
revised.dfs <- subset(revised.dfs, subchallenge == subchallenge.run)
revised.dfs <- revised.dfs[, !(colnames(revised.dfs) == "subchallenge.run")]

orig.nrows <- nrow(res.all)
res.all <- merge(as.data.frame(res.all), as.data.frame(revised.dfs), all.x = TRUE)
new.nrows <- nrow(res.all)
if(orig.nrows != new.nrows) { stop(paste0("Dimension changed from ", orig.nrows, " rows to ", new.nrows, " rows\n")) }
flag <- !is.na(res.all$revised.orig.ratio)
if(any(flag)) {
  res.all[flag, "prediction"] <- res.all[flag, "prediction"] * res.all[flag, "revised.orig.ratio"]
}

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
## REGGEN_LAB is so bad it is skewing the results; leave it off.

exclude.methods <- c("REGGEN_LAB")
exclude.methods <- c("")
## Nah. Instead use annotations to exclude score-based approaches
if(FALSE) {
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
}

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

method.anno <- get.method.annotations()

my.zero.format <- function(x) {
    ifelse(x == 0, 0, x)
}


do.sample.level.analysis <-
    function(res.input, method.anno.round, dataset.annotation,
             method.name.col, 
             model.id.col, subchallenge.col, measured.col, cell.type.col,
             dataset.name.col, sample.id.col, prediction.col,
             round.col, round = "latest", y.sz = 16,
             postfix) {

        submitter.tbl <- unique(res.input[, c(method.name.col, round.col, subchallenge.col), drop = FALSE])
        or <- order(submitter.tbl[, round.col])
        submitter.tbl <- submitter.tbl[or, ]
        submitter.tbl[, round.col] <- as.character(submitter.tbl[, round.col])
        flag <- submitter.tbl[, round.col] %in% c("latest", as.character(round))
        submitter.tbl <- submitter.tbl[flag, ]
        flag <- !duplicated(submitter.tbl[, c(method.name.col, subchallenge.col)], fromLast = FALSE)
        submitter.tbl <- submitter.tbl[flag, ]

        sub.title <- "NA"
        if(round == "1") {
            sub.title <- "First Submission"
        } else if(round == "2") {
            sub.title <- "Up To Second Submission"
        } else if(round == "3") {
            sub.title <- "Up To Third Submission"        
        } else if(round == "latest") {
            sub.title <- "Up To Final Submission"                
        }
        

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
							 ## Nah -- filter this below.
                                                         ## if(df[1, method.name.col] %in% names(deconv.fraction.methods)) {
                                                             val <- sqrt(mean((pred - obs)^2))
                                                         ## } 
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

        ## Plot beeswarm of individuals
        g.swarm.no.color <-
            llply(sub.challenges,
                  .fun = function(sub.challenge) {
                                tbl <- metric.all.res[[sub.challenge]]
                                m.tbl <- melt(tbl)
                                colnames(m.tbl) <- c(method.name.col, dataset.name.col, sample.id.col, "metric", "val")
                                ## m.tbl <- na.omit(m.tbl)
                                m.tbl <- merge(m.tbl, dataset.annotation)
				cat(paste0(sub.challenge, ": # unique samples = ", length(unique(m.tbl[, sample.id.col])), "\n"))
                                decreasing <- FALSE
				metric.order <- "Pearson"
                                order.by.fun <- mean
                                order.by.fun <- median	
                                sm <- ddply(subset(m.tbl, metric == metric.order), .variables = c(method.name.col),
				            .fun = function(df) {
						     data.frame("val" = order.by.fun(df[, "val"]))
						   })

                                flag <- is.na(method.anno.round[, subchallenge.col]) | (as.character(method.anno.round[, subchallenge.col]) == sub.challenge)
				method.anno.round.sc <- method.anno.round[flag, ]
				## Exclude score-based annotations
				method.anno.round.sc <- subset(method.anno.round.sc, Output != "Score")
				## Exclude some methods
				flag <- method.anno.round.sc[, method.name.col] %in% exclude.methods
				method.anno.round.sc <- method.anno.round.sc[!flag, ]
				sm[, method.name.col] <- as.character(sm[, method.name.col])
				for(col in colnames(method.anno.round.sc)) { method.anno.round.sc[, col] <- as.character(method.anno.round.sc[, col]) }

                                m.tbl <- merge(m.tbl, method.anno.round.sc, by = method.name.col, all = FALSE)
                                tbl <- merge(tbl, method.anno.round.sc, by = method.name.col, all = FALSE)
                                sm <- merge(sm, method.anno.round.sc, by = method.name.col, all = FALSE)				
                                o <- order(sm[, "val"], decreasing = decreasing)
                                lvls <- sm[o, method.name.col]
                                m.tbl[, method.name.col] <- factor(m.tbl[, method.name.col], levels = lvls)
                                tbl[, method.name.col] <- factor(tbl[, method.name.col], levels = lvls)
				sm[, method.name.col] <- factor(sm[, method.name.col], levels = lvls)				
				if(FALSE) {
                                  g <- ggplot()
				
                                  g <- g + geom_boxplot(data = m.tbl, aes_string(x = method.name.col, y = "val"))
                                  ## g <- g + facet_wrap("metric", scales = "free_y")
  				  g <- g + facet_wrap(~metric, scales = "free_x")
				  g <- g + coord_flip()
                                  ## g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
                                  g <- g + ylab("") + xlab("")
				}

                                sz <- 16
				title.sz <- 18
                                g1 <- ggplot(data = tbl)
                                g1 <- g1 + geom_boxplot(aes_string(x = method.name.col, y = "Pearson"))
                                g1 <- g1 + coord_flip()
                                g1 <- g1 + xlab("Method")
                                g1 <- g1 + ylab("Pearson Correlation")
                                g1 <- g1 + theme(text = element_text(size=sz), axis.text.y = element_text(size=y.sz), title = element_text(size = title.sz))
                                g1 <- g1 + theme(axis.title.y = element_blank())
                                g1 <- g1 + scale_y_continuous(labels = my.zero.format, limits = c(-1, 1))
				
                                g2 <- ggplot(data = tbl)
                                g2 <- g2 + geom_boxplot(aes_string(x = method.name.col, y = "Spearman"))
                                g2 <- g2 + coord_flip()
                                g2 <- g2 + xlab("Method")
                                g2 <- g2 + ylab("Spearman Correlation")
                                g2 <- g2 + theme(text = element_text(size=sz))    
                                g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                                 axis.ticks.y = element_blank())
                                g2 <- g2 + scale_y_continuous(labels = my.zero.format, limits = c(-1, 1))

                                g3 <- ggplot(data = tbl)
                                g3 <- g3 + geom_boxplot(aes_string(x = method.name.col, y = "RMSE"))
                                g3 <- g3 + coord_flip()
                                g3 <- g3 + xlab("Method")
                                g3 <- g3 + ylab("RMSE")
                                g3 <- g3 + theme(text = element_text(size=sz))    
                                g3 <- g3 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                                 axis.ticks.y = element_blank())
				if(sub.challenge == "fine") {
				  g3 <- g3 + scale_y_continuous(labels = my.zero.format, limits = c(0, 0.3))
				} else {
				  g3 <- g3 + scale_y_continuous(labels = my.zero.format, limits = c(0, 0.3))
				}


                                tmp <- sm[, c(method.name.col, "Output", "Method")]
                                ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))

                                full.plot <- ret[["full.plot"]]
                                for.first.legend <- ret[["legends"]][["Method"]]
                                for.second.legend <- ret[["legends"]][["Output"]]

                                leg1.just <- 1
				## if(sub.challenge == "coarse") { leg1.just <- 0.95 }
                                leg1 <- get_legend(for.first.legend + theme(legend.justification=c(0,leg1.just)))
				leg2 <- get_legend(for.second.legend + theme(legend.justification=c(0,0.7)))
                                ## legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
                                title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge (", sub.title, ")")    				

                                ## plot_row <- plot_grid(g, full.plot, legs, nrow=1, align="h", axis = "b", rel_widths = c(3, 0.5, 0.5))
				## plot_row <- plot_grid(g1, g2, g3, full.plot, legs, nrow=1, align="h", axis = "b", rel_widths = c(5, 3, 3, 0.75, 0.75))
				## plot_row <- plot_grid(g1, g2, g3, full.plot, leg1, leg2, nrow=1, align="h", axis = "b", rel_widths = c(5, 3, 3, 0.75, 0.75, 0.75))
				plot_row_tmp <- plot_grid(g1, g2, g3, full.plot, nrow = 1, align="h", axis = "b", rel_widths = c(5, 3, 3, 0.75))
				plot_row <- plot_grid(plot_row_tmp, leg1, leg2, nrow = 1, rel_widths = c(11.75, 0.75, 0.75))
                                g.with.legs <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

                                g.with.legs
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
                         "g.summaries" = g.summaries, "g.swarm" = g.swarm, "g.swarm.no.color" = g.swarm.no.color)
        return(ret.list)

}

rounds <- c("1", "2", "3")
names(rounds) <- rounds

results <- list()
for(round in rounds) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    method.anno.round <- get.round.specific.annotations(method.anno, round)

    if(FALSE) {
        res.input <- res.all
        round.col <- "submission"
        round <- "2"
        postfix <- paste0("-round-", round)
    }
    
    results[[round]] <- do.sample.level.analysis(res.all, method.anno.round, dataset.annotation,
                                                 method.name.col,
                                                 model.id.col, subchallenge.col, measured.col, cell.type.col,
                                                 dataset.name.col, sample.id.col, prediction.col,
                                                 round.col = "submission", round = round, y.sz = 10,
                                                 postfix)

    sub.title <- "NA"
    if(round == "1") {
        sub.title <- "First Submission"
    } else if(round == "2") {
        sub.title <- "Up To Second Submission"
    } else if(round == "3") {
        sub.title <- "Up To Third Submission"        
    } else if(round == "latest") {
        sub.title <- "Up To Final Submission"                
    }
        
    g.sum.coarse <- results[[round]][["g.summaries"]][["coarse"]]
    title <- paste0("Coarse-Grained (", sub.title, ")")    
    g.sum.coarse <- g.sum.coarse + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    g.sum.fine <- results[[round]][["g.summaries"]][["fine"]]
    title <- paste0("Fine-Grained (", sub.title, ")")        
    g.sum.fine <- g.sum.fine + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    png(paste0(figs.dir, "/sample-level-metric-summary", postfix, ".png"))
    g <- plot_grid(g.sum.coarse, g.sum.fine, nrow = 2, labels = "AUTO") 
    print(g)
    ## grid.arrange(g.sum.coarse, g.sum.fine)
    d <- dev.off()

    g.swarm.coarse <- results[[round]][["g.swarm.no.color"]][["coarse"]]
    title <- paste0("Coarse-Grained (", sub.title, ")")    
##    g.swarm.coarse <- g.swarm.coarse + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    g.swarm.fine <- results[[round]][["g.swarm.no.color"]][["fine"]]
    title <- paste0("Fine-Grained (", sub.title, ")")        
##    g.swarm.fine <- g.swarm.fine + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    png(paste0(figs.dir, "/sample-level-metric-swarm", postfix, ".png"), width = 2 * 480)
    g <- plot_grid(g.swarm.coarse, g.swarm.fine, nrow = 2, labels = "AUTO") 
    print(g)
    ## grid.arrange(g.swarm.coarse, g.swarm.fine)
    d <- dev.off()

}

## Calculate difference relative to top-performer across each submission and each metric

## Exclude these two outliers, which are throwing off the stats
exclude.from.rmse <- c("Patrick", "NYIT_glomerular")

res.tbl <-
  ldply(rounds,
        .fun = function(round) {
	         tbl.sc <-
		   ldply(sub.challenges,
		         .fun = function(sc) {
                                  tbl <- results[[round]][["metric.all.res"]][[sc]]
                                  m.tbl <- melt(tbl)
                                  colnames(m.tbl) <- c(method.name.col, dataset.name.col, sample.id.col, "metric", "val")

                                  method.anno.round <- get.round.specific.annotations(method.anno, round)

                                  flag <- is.na(method.anno.round[, subchallenge.col]) | (as.character(method.anno.round[, subchallenge.col]) == sc)
    				  method.anno.round.sc <- method.anno.round[flag, ]
  				  ## Exclude score-based annotations
				  method.anno.round.sc <- subset(method.anno.round.sc, Output != "Score")
				  ## Exclude some methods
				  flag <- method.anno.round.sc[, method.name.col] %in% exclude.methods
				  method.anno.round.sc <- method.anno.round.sc[!flag, ]
                                  method.anno.round.sc[, method.name.col] <- as.character(method.anno.round.sc[, method.name.col])
                                  m.tbl <- merge(m.tbl, method.anno.round.sc[, c(method.name.col, "Output", "Method")], by = c(method.name.col))

	                          metrics <- as.character(unique(m.tbl$metric))
	                          names(metrics) <- metrics
                                  tbl.met <-
				    ldply(metrics,
	                                  .fun = function(metric.order) {
                                                   ##postfix <- paste0("-", sc, "-round-", round, "-", metric.order)
						   if(metric.order == "RMSE") {
						     flag <- m.tbl[, method.name.col] %in% exclude.from.rmse
						     m.tbl <- m.tbl[!flag, ]
						   }
	                                           m.met <- subset(m.tbl, metric == metric.order)
	                                           m.met <- na.omit(m.met)
                                                   decreasing <- TRUE
                                                   alt <- "greater"
                                                   if(metric.order == "RMSE") { decreasing <- FALSE; alt <- "less" }
                                                   order.by.fun <- median	
                                                   sm <- ddply(m.met, .variables = c(method.name.col),
                                                               .fun = function(df) {
				                                                     data.frame("val" = order.by.fun(df[, "val"]))
                                                                                   })
                                                   o <- order(sm[, "val"], decreasing = decreasing)
                                                   lvls <- as.character(sm[o, method.name.col])
						   do.anova <- TRUE
						   if(do.anova) {
						     m.met[, method.name.col] <- factor(m.met[, method.name.col], levels = lvls)
						     lm.fit <- lm(val ~ method.name, data = m.met)
                                                     sm <- summary(lm.fit)
                                                     cf <- coef(sm)
                                                     flag <- grepl(rownames(cf), pattern="method")
                                                     ret.df <- cf[flag,]
						     val.med <- ddply(m.met, .variables = c(method.name.col), .fun = function(df) data.frame(val = median(df$val)))
						     colnames(val.med)[1] <- method.name.col
						     rownames(val.med) <- val.med[, method.name.col]
						     val.med <- val.med[lvls, ]
                                                     ret.df <- cbind(variable = rownames(ret.df), val.med = val.med[2:nrow(val.med),"val"], ret.df)
                                                     pval <- pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)
                                                     ret.df <- rbind(ret.df, c("F-statistic", NA, as.numeric(sm$fstatistic[1]), NA, NA, pval))
						     ret.df <- cbind(best = lvls[1], best.med = val.med[1,"val"], ret.df)
						     row.names(ret.df) <- NULL
						   } else {
                                                     names(lvls) <- lvls
                                                     flag <- m.met[, method.name.col] == lvls[1]
                                                     res1 <- m.met[flag, ]
                                                     ret.df <-
                                                       ldply(lvls[2:length(lvls)],
                                                             .fun = function(meth) {
                                                                      flag <- m.met[, method.name.col] == meth
                                                                      res2 <- m.met[flag, ]
                                                                      mer <- merge(res1, res2, by = c(sample.id.col, dataset.name.col), suffixes = c(".x", ".y"))
                                                                      print(head(mer))
                                                                      wt <- wilcox.test(mer$val.x, mer$val.y, alternative = alt, paired = TRUE)
                                                                      data.frame(best = lvls[1], p = wt$p.value)
                                                                    })
                                                     colnames(ret.df)[1] <- "comp"			     
                                                   }
						   ret.df
						 })
                                  colnames(tbl.met)[1] <- "metric"
				  tbl.met
                                })
		 colnames(tbl.sc)[1] <- "sub.challenge"
		 tbl.sc
               })	
colnames(res.tbl)[1] <- "round"

flag <- colnames(res.tbl) == "Pr(>|t|)"
colnames(res.tbl)[flag] <- "p.val"

file <- paste0(figs.dir, "/sample-level-comparison.tsv")
write.table(file = file, res.tbl, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat("Exiting successfully\n")

## Summarize based on ties
p.val.cutoff <- 0.05

res.tbl.sum <-
  ddply(res.tbl,
        .variables = c("round", "sub.challenge", "metric"),
	.fun = function(df) {
	         df$p.val <- as.numeric(df$p.val)
                 ties <- as.character(subset(df, p.val > p.val.cutoff)[, "variable"])
		 ties <- gsub(x=ties, pattern="method.name", replacement="")
		 ties <- sort(ties)
		 ties <- paste0(ties, collapse = ", ")
		 data.frame(best = df[1, "best"], ties = ties)
               })

file <- paste0(figs.dir, "/sample-level-comparison-ties.tsv")
write.table(file = file, res.tbl.sum, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat("Exiting successfully\n")
