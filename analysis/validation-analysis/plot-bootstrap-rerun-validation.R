suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(cowplot))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))
suppressPackageStartupMessages(p_load("reshape2"))

suppressPackageStartupMessages(p_load("xlsx"))
suppressPackageStartupMessages(p_load("patchwork"))
suppressPackageStartupMessages(p_load("data.table"))

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

## Read in the bootstrap-validation-results.rds file
synId <- "syn22951683"
obj <- synGet(synId, downloadFile=TRUE)
results <- readRDS(obj$path)

# Condense all of the bootstrap correlations as Source Data for
# Figures 2 and 3 and Supp Figures S5-S15
all.correlations <-
  ldply(results[c("1","2","3")],
        .fun = function(results.submission) {
                 results.subchallenge <- ldply(results.submission[["bootstrapped.cors"]])
                 colnames(results.subchallenge)[1] <- "sub.challenge"
                 results.subchallenge
               })
colnames(all.correlations)[1] <- "submission"
colnames(all.correlations)[colnames(all.correlations) == "boot.i"] <- "bootstrap.index"
fwrite(all.correlations, file="source-data-figs-2-3.csv.gz", quote=FALSE, sep=",")

# Confirm that we can read the results back in
tmp <- fread("source-data-figs-2-3.csv.gz")

# We don't need these files. We just wanted to output the data in a convenient format.
rm(tmp)
rm(all.correlations)
gc()

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

## summary.fun <- mean
summary.fun <- median


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


## for(round in c("1", "2", "3", "latest")) {
plots <- list()
rounds <- c("latest", "1", "2", "3")
# rounds <- c("1")

## Get method metadata
method.anno <- get.method.annotations()

all.median.bootstrapped.scores <- list()

plot.spearman.distribution <- TRUE

source("validation-analysis-utils.R")

for(round in rounds) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    res.round <- results[[round]][["res.round"]]
    flag <- res.round[, "distribution.type"] == "Random"
    res.round[flag, "distribution.type"] <- "Unconstrained"
    bootstrapped.cors <- results[[round]][["bootstrapped.cors"]]
    bootstrapped.scores <- results[[round]][["bootstrapped.scores"]]
    mean.bootstrapped.scores <- results[[round]][["mean.bootstrapped.scores"]]

    ## Recalculate bootstrapped score summary using median
    # We only use median in the boxplot of pearson (implicitly)
    # and in the barplot of spearman (explicitly and for consistency with the pearson boxplot)

    median.bootstrapped.scores <-
        llply(bootstrapped.scores, .parallel = TRUE,
              .fun = function(df) {
                  ## Average over bootstraps
                  df <- ddply(df,
                              .variables = c(method.name.col),
                              .fun = function(tmp) {
                                  data.frame(pearson = summary.fun(tmp$pearson), spearman = summary.fun(tmp$spearman), rmse = summary.fun(tmp$rmse), pearson.fc = summary.fun(tmp$pearson.fc))
                              })
                  o <- order(df$pearson, decreasing = TRUE)
                  df <- df[o,]
                  df
              })
    

    all.median.bootstrapped.scores[[round]] <- median.bootstrapped.scores

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

    ## Need to take annotation for latest round if there isn't one for current round
    method.anno.round <- get.round.specific.annotations(method.anno, round)
	   
    print(head(method.anno.round))
    
    plots[[round]] <- plot.bootstrap.analysis(res.round, bootstrapped.scores, mean.bootstrapped.scores, median.bootstrapped.scores,
                                              means.by.cell.type.method,
                                              means.over.dataset, method.anno.round,
                                              postfix, plot.spearman.distribution = plot.spearman.distribution)
    
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

saveRDS(all.median.bootstrapped.scores, "median-bootstrapped-scores.rds")

for(round in c("1")) {
  tbl <- plots[[round]][["merged.all.means"]]
  tbl <-
    ddply(tbl, .variables = c("cell.type"),
          .fun = function(df) {
	           df <- df[order(df$cor, decreasing=TRUE),]
		   df$diff <- 0
		   df[1:(nrow(df)-1), "diff"] <- df[1:(nrow(df)-1),"cor"] - df[2:(nrow(df)),"cor"]
		   df$comparator.corr <- 0
		   df[1:(nrow(df)-1), "comparator.corr"] <- df[2:(nrow(df)),"cor"]
		   df$comparator <- ""
		   df[1:(nrow(df)-1), "comparator"] <- df[2:(nrow(df)),method.name.col]		   		   
		   df[!is.na(df$cor) & df$cor == max(df$cor, na.rm=TRUE),,drop=F]
		 })
  tbl <- tbl[order(tbl$diff, decreasing=TRUE),]

  cat(paste0("Round ", round, ": ", length(unique(tbl[, method.name.col])), " were top performer in at least one of the ", length(unique(tbl[, cell.type.col])), "\n"))

  file <- paste0(figs.dir, "rerun-validation-round-", round, "-top-method-per-cell-type.tsv")
  write.table(file = file, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}

comparators <- get.comparators.cap()

for(round in c("1")) {
  tbl <- plots[[round]][["merged.all.means"]]
  tbl <-
    ddply(tbl, .variables = c("cell.type"),
          .fun = function(df) {
	           df <- df[order(df$cor, decreasing=TRUE),]
		   flag <- ( !is.na(df$cor) & df$cor == max(df$cor, na.rm=TRUE) ) | ( df[, method.name.col] %in% comparators )
		   df <- df[flag,,drop=F]
		   df$diff <- 0
		   df[1:(nrow(df)-1), "diff"] <- df[1:(nrow(df)-1),"cor"] - df[2:(nrow(df)),"cor"]
		   df$comparator.corr <- 0
		   df[1:(nrow(df)-1), "comparator.corr"] <- df[2:(nrow(df)),"cor"]
		   df$comparator <- ""
		   df[1:(nrow(df)-1), "comparator"] <- df[2:(nrow(df)),method.name.col]		   		   
		   print(df[1,"cell.type"])
		   print(df)
		   df[!is.na(df$cor) & df$cor == max(df$cor, na.rm=TRUE),,drop=F]
		 })
  tbl <- tbl[order(tbl$diff, decreasing=TRUE),]

  cat(paste0("Round ", round, ": ", length(unique(tbl[, method.name.col])), " were top performer in at least one of the ", length(unique(tbl[, cell.type.col])), "\n"))

  file <- paste0(figs.dir, "rerun-validation-round-", round, "-top-method-per-cell-type-rel-comparator.tsv")
  write.table(file = file, tbl, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
}

for(round in c("1")) {
    for(sub.challenge in sub.challenges) {
        # Note that this is using the _mean_ of bootstrapped scores in results[[round]],
        # not the _median_ of bootstrapped scores in plots[[round]].
	# This is for consistency with other plots (e.g., the heatmaps, where we
	# report the means). We only use median in the boxplot of pearson (implicitly)
	# and in the barplot of spearman (explicitly and for consistency with the pearson boxplot)
        # Actually, here we use the mean in plots[[round]] as this has an Output column
        # appended to the mean in results[[round]]
        # tbl <- plots[[round]][["median.bootstrapped.scores"]][[sub.challenge]]
        # tbl <- results[[round]][["mean.bootstrapped.scores"]][[sub.challenge]]
        tbl <- plots[[round]][["mean.bootstrapped.scores"]][[sub.challenge]]
        # Exclude the ensemble method
        flag <- tbl[, method.name.col] %in% c("consensus rank", "ensemble")
        tbl <- tbl[!flag, ]
	tbl <- as.data.frame(table(na.omit(tbl$Method)))
	colnames(tbl) <- c("Method", "Freq")
	o <- order(tbl$Freq)
	tbl <- tbl[o, ]
	cat(paste0("Round ", round, " ", sub.challenge, " method count:\n"))
	print(tbl)

        tbl <- plots[[round]][["mean.bootstrapped.scores"]][[sub.challenge]]
        # Exclude the ensemble method
        flag <- tbl[, method.name.col] %in%  c("consensus rank", "ensemble")
        cat(paste0("Excluded ensemble method: ", any(flag), "\n"))
        tbl <- tbl[!flag, ]
	mean.tbl <- ddply(tbl, .variables = c("Output"),
	             .fun = function(df) {
		              data.frame(pearson = mean(df$pearson, na.rm=TRUE),
			                 spearman = mean(df$spearman, na.rm=TRUE))
			    })
	colnames(mean.tbl) <- c("Output", "pearson", "spearman")
	o <- order(mean.tbl$pearson)
	mean.tbl <- mean.tbl[o, ]
	cat(paste0("Round ", round, " ", sub.challenge, " method count:\n"))
	print(mean.tbl)
	pwt <- pairwise.wilcox.test(tbl$pearson, tbl$Output, p.adjust.method = "BH")
	print(pwt)

	
    }
}

dataset.scores <- list()
non.nas <- list()
for(round in rounds) {
    dataset.scores[[round]] <- list()
    non.nas[[round]] <- list()
    for(sub.challenge in sub.challenges) {
        df <- results[[round]][["means.over.bootstrap"]][[sub.challenge]][["pearson"]]
        methods.with.nas <- unique(subset(df, is.na(cor))[, method.name.col])
        flag <- !(df[, method.name.col] %in% methods.with.nas)
        means <- ddply(df[flag, ], .variables = c(method.name.col, dataset.name.col),
                       .fun = function(sub) data.frame(cor = mean(sub$cor, na.rm=FALSE)))
        anno <- unique(results[[round]][["res.round"]][, c(dataset.name.col, "mixture.type", "distribution.type")])
        means <- merge(means, anno)
        non.nas[[round]][[sub.challenge]] <- merge(df[flag, ], anno)
        dataset.scores[[round]][[sub.challenge]] <- means
    }
}

# rounds <- list("1" = "1", "2" = "2", "3" = "3")
names(rounds) <- rounds
metrics <- list("pearson" = "pearson", "spearman" = "spearman", "rmse" = "rmse", "pearson.fc" = "pearson.fc")
lm.fits <- ldply(rounds,
             .fun = function(round) {
                 ret2 <- ldply(sub.challenges,
                               .fun = function(sub.challenge) {
                                   ret1 <- ldply(metrics,
                                                 .fun = function(metric) {
                                                     df <- results[[round]][["means.over.bootstrap"]][[sub.challenge]][[metric]]
                                                     methods.with.nas <- unique(subset(df, is.na(cor))[, method.name.col])
                                                     flag <- !(df[, method.name.col] %in% methods.with.nas)
                                                     df <- df[flag, ]
                                                     anno <- unique(results[[round]][["res.round"]][, c(dataset.name.col, "mixture.type", "distribution.type")])
                                                     df <- merge(df, anno)
                                                     lm.fit <- lm(cor ~ mixture.type + distribution.type + method.name + cell.type, data = df)
						     sm <- summary(lm.fit)
                                                     cf <- coef(sm)
                                                     flag <- grepl(rownames(cf), pattern="distribution") | grepl(rownames(cf), pattern="mixture")
                                                     ret.df <- cf[flag,]
                                                     ret.df <- cbind(variable = rownames(ret.df), ret.df)
						     pval <- pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)
						     ret.df <- rbind(ret.df, c("F-statistic", as.numeric(sm$fstatistic[1]), NA, NA, pval))
                                                 })
                                   colnames(ret1)[1] <- "metric"
                                   ret1
                               })
                 colnames(ret2)[1] <- subchallenge.col
                 ret3 <- ldply(metrics,
                               .fun = function(metric) {
                                   df <- rbind(results[[round]][["means.over.bootstrap"]][["coarse"]][[metric]],
                                               results[[round]][["means.over.bootstrap"]][["fine"]][[metric]])
                                   methods.with.nas <- unique(subset(df, is.na(cor))[, method.name.col])
                                   flag <- !(df[, method.name.col] %in% methods.with.nas)
                                   df <- df[flag, ]
                                   df <- ddply(df,
                                               .variables = c(method.name.col, cell.type.col, dataset.name.col),
                                               .fun = function(sub) data.frame(cor = mean(sub$cor)))
                                   anno <- unique(results[[round]][["res.round"]][, c(dataset.name.col, "mixture.type", "distribution.type")])
                                   df <- merge(df, anno)
                                   lm.fit <- lm(cor ~ mixture.type + distribution.type + method.name + cell.type, data = df)
				   sm <- summary(lm.fit)
                                   cf <- coef(sm)
                                   flag <- grepl(rownames(cf), pattern="distribution") | grepl(rownames(cf), pattern="mixture")
                                   ret.df <- cf[flag,]
                                   ret.df <- cbind(variable = rownames(ret.df), ret.df)
		                   pval <- pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3],lower.tail=FALSE)
		                   ret.df <- rbind(ret.df, c("F-statistic", as.numeric(sm$fstatistic[1]), NA, NA, pval))
                               })
                 colnames(ret3)[1] <- "metric"
                 ret3 <- cbind("merged", ret3)
                 colnames(ret3)[1] <- subchallenge.col
                 rbind(ret2, ret3)
             })
colnames(lm.fits)[1] <- round.col
sig.cutoff <- 0.01
lm.fits[,8] <- as.numeric(as.character(lm.fits[,8]))
flag <- lm.fits[, 8] < sig.cutoff

plot.mixture.distribution.effect <- function(df) {
    flag <- df[, "distribution.type"] == "Random"
    df[flag, "distribution.type"] <- "Unconstrained"
    g <- ggplot(data = df, aes_string(x = dataset.name.col, y = "cor", 
                                                colour = "mixture.type",
                                                linetype = "distribution.type"))
    g <- g + geom_boxplot()
    ## , scales = "free_y")
    g <- g + facet_wrap(method.name.col)
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g <- g + scale_linetype_discrete(name = "Distribution")
    g <- g + scale_colour_discrete(name = "Mixture")
    g <- g + ylab("") + xlab("")
    g 
}


plot.merged.mixture.distribution.effiect <- function(results, round = "1", metric = "pearson") {
    df <- rbind(results[[round]][["means.over.bootstrap"]][["coarse"]][[metric]],
                results[[round]][["means.over.bootstrap"]][["fine"]][[metric]])

    methods.with.nas <- unique(subset(df, is.na(cor))[, method.name.col])
    flag <- !(df[, method.name.col] %in% methods.with.nas)
    df <- df[flag, ]
    df <- ddply(df,
                .variables = c(method.name.col, cell.type.col, dataset.name.col),
                .fun = function(sub) data.frame(cor = mean(sub$cor)))
    anno <- unique(results[[round]][["res.round"]][, c(dataset.name.col, "mixture.type", "distribution.type")])
    df <- merge(df, anno)
    g <- plot.mixture.distribution.effect(df)
    g
}

g <- plot.merged.mixture.distribution.effiect(results, round = "1", metric = "pearson")
g <- g + ylab("Pearson Correlation")
g <- g + ggtitle("Merged Coarse- and Fine-Grained (First Submission)")
g <- g + theme(plot.title = element_text(hjust = 0.5))

file.suffix <- paste0(figs.dir, "fig-validation-round-1-merged-mixture-distribution-effect")
output.plot(g, file.suffix, plot.types = c("pdf", "png"),
            pdf.delta.width = 1, pdf.delta.height = 1,
            png.delta.width = 1, png.delta.height = 1)

#png(paste0(figs.dir, "fig-validation-round-1-merged-mixture-distribution-effect.png"))
#print(g)
#d <- dev.off()

#pdf(paste0(figs.dir, "fig-validation-round-1-merged-mixture-distribution-effect.pdf"))
#print(g)
#d <- dev.off()

file <- paste0(figs.dir, "rerun-validation-mixture-and-distribution-effects.tsv")
write.table(file = file, lm.fits, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

plot.scores.over.rounds <- function(df, comparator.methods) {
    order.round <- "1"
    df.order <- subset(df, Round == order.round)
    o <- order(df.order$pearson)
    lvls <- df.order[o, method.name.col]
    df <- subset(df, Round != "latest")
    lvls <- lvls[lvls %in% df[, method.name.col]]
    bold.labels <- ifelse(lvls %in% comparator.methods, yes = "bold", no = "plain")
print("plot.scores.over.rounds")
print(comparator.methods)
print(lvls)
print(bold.labels)
    df[, method.name.col] <- factor(df[, method.name.col], levels = lvls)
    df$Submission <- df$Round
    g <- ggplot()
    g <- g + geom_point(data = df, aes_string(x = "pearson", y = method.name.col, colour = "Submission"))
    g <- g + xlab("Pearson Correlation") + ylab("Method")
    g <- g + theme(axis.text.y = element_text(face = bold.labels))
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
    llply(sub.challenges, 
          .fun = function(sc) {
                   comparator.methods <- llply(results, .fun = function(df) unique(subset(df[["res.round"]], subchallenge == sc & comparator == TRUE)[, method.name.col]))
                   comparator.methods <- unique(as.vector(unlist(comparator.methods)))
                   plot.scores.over.rounds(all.means[[sc]], comparator.methods)
                 })

if(FALSE) {
    l_ply(sub.challenges,
          .fun = function(subchallenge) {
              g <- g.score.vs.round[[subchallenge]]
              
              file.suffix <- paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", subchallenge)
              output.plot(g, file.suffix, plot.types = c("pdf", "png"),
                          pdf.delta.width = 1, pdf.delta.height = 1,
                          png.delta.width = 1, png.delta.height = 1)


              #png(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", subchallenge, ".png"))
              #print(g)
              #d <- dev.off()

              #pdf(paste0(figs.dir, "rerun-validation-bootstrap-pearson-vs-round-", subchallenge, ".pdf"))
              #print(g)
              #d <- dev.off()
              
          })
}

shift.limit <- function(val) {
    sign(val) * round(abs(val) + 0.05, digits = 1)
}

source("make-validation-performance-figs.R")

save.image(".Rdata.plot.bootstrap")

cat("Exiting successfully\n")
q(status=0)
