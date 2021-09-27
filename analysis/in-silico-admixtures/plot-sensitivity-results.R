
suppressPackageStartupMessages(library(pacman))

suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(tidyr))
suppressPackageStartupMessages(p_load(grid))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(openxlsx))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(ggpubr))
suppressPackageStartupMessages(p_load(ggbeeswarm))
suppressPackageStartupMessages(p_load(cowplot))

## sensitivity / spike-in analysis

## Get the sensitivity / spike-in results
synLogin()

subchallenge.col <- "subchallenge"
measured.col <- "measured"
cell.type.col <- "cell.type"
dataset.name.col <- "dataset.name"
sample.id.col <- "sample.id"
prediction.col <- "prediction"
method.name.col <- "method.name"
round.col <- "submission"
mixture.col <- "mixture.type"

source("../utils.R")

figs.dir <- "figs"
dir.create(figs.dir, showWarnings = FALSE)

url.base <- "https://github.com/Sage-Bionetworks/Tumor-Deconvolution-Challenge/blob/master/analysis/"
this.script <- "plot-sensitivity-analysis.R"
script_url <- paste0(url.base, "/", "in-silico-admixtures", "/", this.script)

synIds <- list("sensitivity-analysis-data" = "syn22953684",
               "sensitivity-analysis-results" = "syn22953003")

tbls <-
    llply(synIds,
          .fun = function(synId) {
              file <- synGet(synId, downloadFile = TRUE)$path
              readRDS(file)
          })

my.format <- function(x) {
    ifelse(x == 0, 0,
    ifelse(x < 0.01, formatC(x, format="e", digits=1),
    ifelse(x < 1, formatC(x, format="f", digits=1, drop0trailing = TRUE),
           round(x))))
}

my.format2 <- function(x) {
    ifelse(x == 0, 0,
    ifelse(x < 0.01, formatC(x, format="e", digits=1),
    ifelse(x < 1, formatC(x, format="f", digits=2, drop0trailing = TRUE),
           round(x))))
}

summary.plots <- list()
rounds <- c("1", "2", "3", "latest")
for(round in rounds) {
    cat(paste0("Doing round ", round, "\n"))

    res.matrices <- tbls[["sensitivity-analysis-results"]][[round]]$res.matrices

    nms <- names(res.matrices)
    names(nms) <- nms
    summary.plots[[round]] <- 
        llply(nms,
              .fun = function(nm) {
                  res <- res.matrices[[nm]]
                  subplots <-
                      dlply(res, .variables = c("mixture.type"),
                            .fun = function(res.ds) {
                                
                                df <- acast(res.ds, as.formula(paste0(method.name.col, " ~ ", cell.type.col)), value.var = "min.diff.prop", fill = NA)
                                res.ds <- reshape2::melt(as.matrix(df))
                                colnames(res.ds) <- c(method.name.col, cell.type.col, "min.diff.prop")
                                
                                id.var <- method.name.col
                                
                                res.ds$label <- as.numeric(as.character(100 * res.ds$min.diff.prop))
                                
                                ## Cast from long to matrix form to introduce NAs for missing values,
                                ## then melt back to long form
                                
                                g <- plot.cell.type.correlation.heatmap(res.ds, show.corr.text = TRUE,
                                                                        id.var = method.name.col, cell.type.var = cell.type.col,
                                                                        cor.var = "label", formatter = my.format,
                                                                        col.summary.fun = "min", cor.type.label = "LoD (Percent)",
                                                                        limits = c(0, 100),
                                                                        order.decreasing = TRUE)
                                return(g)

                                ## limits = c(log2(min(res$measured[res$measured > 0])), log2(1))
                                limits = c(log2(10^-4), log2(1))
                                my.labeller <- function(x) { my.format(100*(2^x)) }
                                g <- g + scale_fill_gradient2("LoD (Percent)", labels = my.labeller, limits = limits, 
                                                              low = "red", high = "blue", mid = "white", na.value = "black")
                                
                                g
                            })
              })

    ## Output the plots
    for(sub.challenge in names(summary.plots[[round]])) {
        plts <- summary.plots[[round]][[sub.challenge]]
        for(mixture.type in names(plts)) {
            o.file <- paste0(figs.dir, "/", "sensitivity-summary-", sub.challenge, "-round-", round, "-", mixture.type, ".png")
            png(o.file)
            g <- plts[[mixture.type]]
            g <- g + ggtitle(paste0(firstup(sub.challenge), "-Grained Sub-Challenge (Round ", round, ", ", mixture.type, " Admixtures", ")"))
            print(g)
            d <- dev.off()
        }
        g.bio <- plts[["Biological"]]
        g.bio <- g.bio + ggtitle("Biological Admixtures")
        g.rand <- plts[["Random"]]
        g.rand <- g.rand + ggtitle("Random Admixtures")
        o.file <- paste0(figs.dir, "/", "sensitivity-summary-", sub.challenge, "-round-", round, "-both-mixture-types.png")
        png(o.file, width = 2 * 480, height = 1 * 480)
        title <- paste0(firstup(sub.challenge), "-Grained Sub-Challenge (Round ", round, ")")
        g <- grid.arrange(g.bio, g.rand, nrow = 1, top = textGrob(title, gp=gpar(fontsize=25)))
        grid.draw(g)
        d <- dev.off()
    }
}

plot.spikein.predictions <- function(df.ds) {
    ct <- as.character(df.ds[1, cell.type.col])
    df.ds <- df.ds[order(df.ds[, measured.col]),]
    df.ds$val <- as.numeric(as.character(df.ds[, measured.col])) * 100
    ## sci <- formatC(as.numeric(as.character(df.ds$measured))*100, format="e", digits=2)
    sci <- my.format2(df.ds$val)
    df.ds$val <- factor(sci, levels = unique(sci))
    cmps <- as.data.frame(compare_means(prediction ~ val,  data = df.ds))
    ## cmps <- cmps[cmps$group1 == "0.00e+00",]
    cmps <- cmps[cmps$group1 == "0",]
    display.p.cutoff <- 0.1
    display.p.cutoff <- 0.01
    cmps <- cmps[cmps$p < display.p.cutoff,]
    my_comparisons <-
        llply(1:nrow(cmps),
              .fun = function(i) c(cmps[i,2], cmps[i,3]))
                                  
    g <- ggplot(data = df.ds, aes_string(x = "val", y = prediction.col))
    g <- g + geom_boxplot()
    ## Remove label = ..p.signif.. to plot the actual pvalues
    g <- g + stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, size = 5, vjust = 0.5)
    g <- g + geom_beeswarm()
    g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    xlab <- paste0(fix.string(ct), " Spike-in (%)")
    g <- g + xlab(xlab) + ylab("Score")
    g
}

expts.to.plot <-
    list(
        c("EPIC", "B", "coarse", "Random"),
        c("CIBERSORTx", "CD4 T", "coarse", "Random"),
        c("Aginome-XMU", "CD4 T", "coarse", "Random"),
        c("DA_505", "memory CD4 T", "fine", "Random")
    )
expts.to.plot <- t(as.data.frame(expts.to.plot))
rownames(expts.to.plot) <- NULL
colnames(expts.to.plot) <- c(method.name.col, cell.type.col, subchallenge.col, mixture.col)

spike.in.plots <- list()

for(round in rounds) {
    cat(paste0("Doing round ", round, "\n"))

    res <- tbls[["sensitivity-analysis-results"]][[round]]$res

    res <- merge(res, expts.to.plot)

    spike.in.plots[[round]] <- 
        dlply(res,
              .variables = c(method.name.col, cell.type.col, subchallenge.col, mixture.col),
              .fun = function(df) {
                  expected.num <- 10
                  if(!all(table(df[, c(cell.type.col, "measured")]) == expected.num)) {
                      print(table(df[, c(cell.type.col, "measured")]))
                      stop(paste0("Some spike-in levels did not sum to ", expected.num, "\n"))
                  }
                  
                  mn <- as.character(df[1, method.name.col])
                  ct <- as.character(df[1, cell.type.col])
                  sc <- as.character(df[1, subchallenge.col])
                  mx <- as.character(df[1, mixture.col])
                  g <- plot.spikein.predictions(df)
                  g <- g + ggtitle(paste0(mn, ": ", ct, "\n", firstup(sc), "-Grained Sub-Challenge (", mx, " Admixtures)"))
                  o.file <- paste0(figs.dir, "/", "sensitivity-spikein-", make.names(mn), "-",
                                   make.names(ct), "-", sc, "-", mx, "-round-", round, ".png")
                  png(o.file, width = 2 * 480)
                  print(g)
                  d <- dev.off()
                  return(g)
              })
}


o.file <- paste0(figs.dir, "/", "sensitivity-spikein-and-summary.png")
png(o.file, width = 2 * 480, height = 2 * 480)
g1 <- spike.in.plots[["1"]][["Aginome-XMU.CD4 T.coarse.Random"]]
g2 <- summary.plots[["1"]][["coarse"]][["Random"]]
g3 <- summary.plots[["1"]][["fine"]][["Random"]]
bottom_row <- plot_grid(g2, g3, labels = c("B", "C"))
g <- plot_grid(g1, bottom_row, labels = c("A", ""), ncol =1)
## g <- plot_grid(g1, g2, g3, labels = c("A", "B", "C"))
## g <- grid.arrange(g1, g2, g3)
## grid.draw(g)
print(g)
d <- dev.off()


cat("Exiting successfully\n")
q(status=0)

