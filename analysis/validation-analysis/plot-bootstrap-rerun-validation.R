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

## See https://stackoverflow.com/questions/27803710/ggplot2-divide-legend-into-two-columns-each-with-its-own-title
plot.anno.heatmap.with.multiple.legends <-
    function(df, id.col, anno.columns, anno.pals) {

        suppressPackageStartupMessages(p_load("RColorBrewer"))
        df <- df[, c(id.col, anno.columns)]

        ## Assume annotations are characters
        ## NB: id.col is a factor
        for(col in c(anno.columns)) {
            df[, col] <- as.character(df[, col])
        }

        columns <- 1:length(anno.columns)
        names(columns) <- anno.columns

        color.vecs <-
            llply(columns,
                  .fun = function(idx) {
                      anno.col <- anno.columns[idx]
                      vec <- unique(df[, anno.col])
                      len <- length(vec)
                      colors <- brewer.pal(len, anno.pals[idx])
                      names(colors) <- vec
                      colors
                  })

        all.colors <- Reduce("c", color.vecs)
        names(all.colors) <- Reduce("c", unlist(lapply(color.vecs, names)))

        names(anno.columns) <- anno.columns
        anno.df <- ldply(anno.columns,
                     .fun = function(anno.col) {
                         data.frame(val = df[, anno.col], id = df[, id.col])
                     })
        colnames(anno.df)[1] <- "type"
        print(anno.df)
        full.plot <-
            ggplot(anno.df, aes(y = id, x = type, fill = val)) + geom_tile() +
            scale_fill_manual(values = all.colors) +
            theme(legend.position="none")

        full.plot <- full.plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                       axis.ticks.y = element_blank(), text = element_text(size = 18),
                                       axis.text.x = element_text(angle = 45, hjust = 1),
                                       axis.title.x = element_blank())
        
        
        legends <-
            llply(anno.columns,
                  .fun = function(anno.col) {
                      flag <- anno.df$type == anno.col
                      g <- ggplot(anno.df[flag, ], aes_string(x = "id", y = "type", fill = "val"))
                      g <- g + geom_tile()
                      g <- g + scale_fill_manual(values = all.colors, name = anno.col)
                  })

        return(list("full.plot" = full.plot, "legends" = legends))
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
             means.over.dataset, method.anno.round,
             postfix) {

        top.performers <- c("Aginome-XMU", "DA_505", "mitten_TDC19", "Biogem")
        priority.methods <- unique(c(top.performers, unique(subset(res, comparator==TRUE)[, method.name.col])))

        for(sub.challenge in sub.challenges) {
            print(sub.challenge)
            print(head(method.anno.round))
            print(subchallenge.col)
            flag <- is.na(method.anno.round[, subchallenge.col]) | ( method.anno.round[, subchallenge.col] == sub.challenge )
            print(flag)
            method.anno.round.sc <- method.anno.round[flag, c(method.name.col, "Output", "Method")]
            bootstrapped.scores[[sub.challenge]] <-
                merge(bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
            mean.bootstrapped.scores[[sub.challenge]] <-
                merge(mean.bootstrapped.scores[[sub.challenge]], method.anno.round.sc, all.x = TRUE)
        }
        
        if(FALSE) {
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
                
                tmp <- scores[, c(method.name.col, "Output", "Method")]
                ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))
                
                full.plot <- ret[["full.plot"]]
                for.first.legend <- ret[["legends"]][["Method"]]
                for.second.legend <- ret[["legends"]][["Output"]]
                
                legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
                
                ## pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,1,0.5,0.5))
                
                boxplots[[paste0(sub.challenge,  "-pearson")]] <- g1
                boxplots[[paste0(sub.challenge,  "-spearman")]] <- g2
                boxplots[[paste0(sub.challenge,  "-anno")]] <- full.plot
                boxplots[[paste0(sub.challenge,  "-legend")]] <- legs
                
                title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge")
                round.text <- make.round.text(round)
                title <- paste0(title, " (", round.text, ")")
                ## png(paste0(figs.dir, "rerun-validation-score-boxplots-", sub.challenge, postfix, ".png"), width = 2 * 480)
                ## g <- grid.arrange(g1, g2, nrow=1, top = textGrob(title, gp = gpar(fontsize = 25)))
                ## d <- dev.off()
            }
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
            ## g1 <- g1 + ylab("Pearson Correlation")
            g1 <- g1 + ylab("Pearson")
            g1 <- g1 + ylim(c(-0.25, 1))
            g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20))

            g2 <- ggplot(data = mean.scores)
            g2 <- g2 + geom_col(aes_string(x = method.name.col, y = "spearman"), fill = "#E69F00")
            g2 <- g2 + coord_flip()
            g2 <- g2 + xlab("Method")
            ## g2 <- g2 + ylab("Spearman Correlation")
            g2 <- g2 + ylab("Spearman")
            g2 <- g2 + ylim(c(-0.25, 1))            
            g2 <- g2 + theme(text = element_text(size=18))    
            g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                             axis.ticks.y = element_blank())

            tmp <- scores[, c(method.name.col, "Output", "Method")]            
            ret <- plot.anno.heatmap.with.multiple.legends(tmp, "method.name", c("Method", "Output"), c("Set3", "Set1"))
            
            full.plot <- ret[["full.plot"]]
            for.first.legend <- ret[["legends"]][["Method"]]
            for.second.legend <- ret[["legends"]][["Output"]]
            
            legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(2,1))
            
            ## pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,1,0.5,0.5))
            
            barplots[[paste0(sub.challenge,  "-pearson")]] <- g1
            barplots[[paste0(sub.challenge,  "-spearman")]] <- g2            
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
        method.levels <- list()
        cell.type.levels <- list()
        for(sub.challenge in sub.challenges) {
            means <- means.by.cell.type.method[[sub.challenge]][["pearson"]]

            method.levels[[sub.challenge]] <-
                calculate.method.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")

            cell.type.levels[[sub.challenge]] <-
                calculate.cell.type.levels(means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
            
            
            g <- plot.cell.type.correlation.heatmap(means, show.corr.text = TRUE,
                                                    id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                                    method.levels = method.levels[[sub.challenge]],
                                                    cell.type.levels = cell.type.levels[[sub.challenge]])
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

        method.levels[["merged"]] <-
            calculate.method.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        
        cell.type.levels[["merged"]] <-
            calculate.cell.type.levels(all.means, id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor")
        
        g <- plot.cell.type.correlation.heatmap(all.means, show.corr.text = TRUE,
                                                id.var = method.name.col, cell.type.var = cell.type.col, cor.var = "cor",
                                                method.levels = method.levels[["merged"]],
                                                cell.type.levels = cell.type.levels[["merged"]])
        heatmaps[["merged"]] <- g

        cat("Creating merged means.over.dataset\n")
        coarse.means <- means.over.dataset[["coarse"]][["pearson"]]
        fine.means <- means.over.dataset[["fine"]][["pearson"]]        
        all.means <- rbind(coarse.means, fine.means)

        all.means <-
            ddply(all.means, .variables = c(method.name.col, cell.type.col, "boot.i"),
                  .fun = function(df) {
                      data.frame(cor = mean(df$cor))
                  })

        cat(paste0("Plotting strip plots\n"))
        nms <- list("coarse" = "coarse", "fine" = "fine", "coarse-priority" = "coarse-priority", "fine-priority" = "fine-priority",
                    "merged" = "merged", "merged-priority" = "merged-priority")
        strip.plots <-
            llply(nms,
                  .parallel = TRUE,
                  .fun = function(nm) {
                      sub.challenge <- NA
                      df <- NULL
                      entry <- NULL
                      if(grepl(nm, pattern="coarse")) {
                          entry <- "coarse"
                          df <- means.over.dataset[[entry]][["pearson"]]
                      }
                      if(grepl(nm, pattern="fine")) {
                          entry <- "fine"
                          df <- means.over.dataset[[entry]][["pearson"]]                          
                      }
                      if(grepl(nm, pattern="merged")) {
                          entry <- "merged"
                          df <- all.means
                      }

                      g <- NULL
                      if(grepl(nm, pattern="priority")) {
                          flag <- df[, method.name.col] %in% priority.methods
                          g <- plot.strip.plots(df[flag, ], id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                                method.levels = method.levels[[entry]],
                                                cell.type.levels = cell.type.levels[[entry]],
                                                label = "Pearson Correlation")
                      } else {
                          g <- plot.strip.plots(df, id.var = method.name.col, cell.type.var = cell.type.col, var = "cor",
                                                method.levels = method.levels[[entry]],
                                                cell.type.levels = cell.type.levels[[entry]],
                                                label = "Pearson Correlation")
                      }
                  })

        ## "boxplots" = boxplots,
        ret.list <- list("mean.bootstrapped.scores" = mean.bootstrapped.scores,
                         "barplots" = barplots,
                         "strip.plots" = strip.plots,
                         "heatmaps" = heatmaps)

        ret.list
        
}

## for(round in c("1", "2", "3", "latest")) {
plots <- list()
rounds <- c("2", "1", "3", "latest")
## rounds <- c("1")

## Get method metadata
synId <- "syn23395242"
obj <- synGet(synId, downloadFile = TRUE)
method.anno <- read.xlsx(obj$path, sheetIndex = 1)

method.anno$method.type <- as.character(method.anno$method.type)
method.anno$output.type <- as.character(method.anno$output.type)
method.anno[, round.col] <- as.character(method.anno[, round.col])
flag <- method.anno[, round.col] == "NA"
method.anno[flag, round.col] <- NA
method.anno[, subchallenge.col] <- as.character(method.anno[, subchallenge.col])
flag <- method.anno[, subchallenge.col] == "NA"
method.anno[flag, subchallenge.col] <- NA

method.rename.list <-
    list("NNLS" = "NNLS",
         "summary" = "SUM",
         "other" = "OTH",
         "other regression" = "REG",
         "unknown" = "UNK",
         "SVR" = "SVR",
         "DNN" = "DNN",
         "ensemble" = "ENS",
         "NMF" = "NMF",
         "probabilistic inference" = "PI")
method.rename.df <- data.frame(method.type = names(method.rename.list), Method = as.character(method.rename.list))

output.rename.list <-
    list("fraction" = "Frac",
         "proportion" = "Prop",
         "normalized.score" = "Norm",
         "score" = "Score")
output.rename.df <- data.frame(output.type = names(output.rename.list), Output = as.character(output.rename.list))

method.anno <- merge(method.anno, method.rename.df)
method.anno <- merge(method.anno, output.rename.df)
method.anno$Output <- as.character(method.anno$Output)
method.anno$Method <- as.character(method.anno$Method)

for(round in rounds) {
    postfix <- paste0("-round-", round)
    cat(paste0("Doing round ", round, "\n"))

    res.round <- results[[round]][["res.round"]]
    flag <- res.round[, "distribution.type"] == "Random"
    res.round[flag, "distribution.type"] <- "Unconstrained"
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

    flag <- is.na(method.anno[, round.col]) | (method.anno[, round.col] == round)
    method.anno.round <- method.anno[flag, ]

    print(head(method.anno.round))
    
    plots[[round]] <- plot.bootstrap.analysis(res.round, bootstrapped.scores, mean.bootstrapped.scores,
                                              means.by.cell.type.method,
                                              means.over.dataset, method.anno.round,
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

for(round in c("1")) {
    for(sub.challenge in sub.challenges) {
        tbl <- plots[[round]][["mean.bootstrapped.scores"]][[sub.challenge]]
	tbl <- as.data.frame(table(na.omit(tbl$Method)))
	colnames(tbl) <- c("Method", "Freq")
	o <- order(tbl$Freq)
	tbl <- tbl[o, ]
	cat(paste0("Round ", round, " ", sub.challenge, " method count:\n"))
	print(tbl)

        tbl <- plots[[round]][["mean.bootstrapped.scores"]][[sub.challenge]]
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

rounds <- list("1" = "1", "2" = "2", "3" = "3")
metrics <- list("pearson" = "pearson", "spearman" = "spearman", "rmse" = "rmse")
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

png(paste0(figs.dir, "fig-validation-round-1-merged-mixture-distribution-effect.png"))
print(g)
d <- dev.off()

file <- paste0(figs.dir, "rerun-validation-mixture-and-distribution-effects.tsv")
write.table(file = file, lm.fits, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

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

save.image(".Rdata.plot.bootstrap")

cat("Exiting successfully\n")
q(status=0)
