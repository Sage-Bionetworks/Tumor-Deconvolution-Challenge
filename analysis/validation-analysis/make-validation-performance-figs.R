
## Count number of participants and submissions
## Read in the rerun predicitons (i.e., where the coarse- and fine-grained datasets are the same)
synId <- "syn22320329"
obj <- synGet(synId, downloadFile=TRUE)
res.all <- read.table(obj$path, sep=",", header=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

all.methods <- sort(unique(res.all[, method.name.col]))
annotated.methods <- sort(unique(method.anno[, method.name.col]))

flag <- !(annotated.methods %in% all.methods)
if(any(flag)) {
    cat(paste0("Annotated methods not in results: ", paste0(annotated.methods[flag], collapse=","), "\n"))
}

flag <- !(all.methods %in% annotated.methods)
if(any(flag)) {
    cat(paste0("Resulted methods without annotations: ", paste0(all.methods[flag], collapse=", "), "\n"))
}

##target.meth <- "REGGEN_LAB"
##source("infer-output-type.R")

## sum(subset(d3, sample.id==sample & (submission==1) & (dataset.name==ds) & (subchallenge==sc))$prediction)
                         
submission.tbl <- unique(res.all[, c(method.name.col, round.col, "comparator", subchallenge.col)])

for(sub.challenge in sub.challenges) {
    tbl <- submission.tbl
    flag <- ( tbl$comparator == FALSE ) & ( tbl[, subchallenge.col] == sub.challenge ) &
        !(tbl[, round.col] == "latest")
    tbl <- tbl[flag, ]
    cat(paste0(sub.challenge, " ", length(unique(tbl[, method.name.col])), " teams had ",
               nrow(unique(tbl[, c(method.name.col, round.col)])), "\n"))
}

print(sort(unique(subset(submission.tbl, comparator == TRUE)[, c(method.name.col)])))

compute.chained.empirical.bayes <- function(summarized.bootstrapped.scores, bootstrapped.scores, sub.challenge, metric = "pearson") {
    summarized.scores <- summarized.bootstrapped.scores[[sub.challenge]]
    summarized.scores <- na.omit(summarized.scores[, c(method.name.col, metric)])
    o <- order(summarized.scores[, metric], decreasing = TRUE)
    summarized.scores <- summarized.scores[o, ]
    scores <- bootstrapped.scores[[sub.challenge]]
    indices <- 1:(nrow(summarized.scores)-1)
    names(indices) <- summarized.scores[indices, method.name.col]
    ret <-
        ldply(indices,
              .fun = function(i) {
                  numerator.id <- summarized.scores[i, method.name.col]                  
                  denominator.id <- summarized.scores[i+1, method.name.col]
                  bf <- calculate.empirical.bayes(scores, col.id = method.name.col,
                                                  numerator.id = numerator.id,
                                                  denominator.id = denominator.id,
                                                  sample.id.cols = "boot.i",
                                                  score.col = metric)
                  data.frame(comparator = denominator.id, comparator.rank = i+1, bayes.factor = bf)
              })
    colnames(ret)[1] <- "method"
    ret
}

chained.empirical.bayes <-
    llply(rounds,
          .fun = function(round) {
              bootstrapped.scores <- results[[round]][["bootstrapped.scores"]]
              # Note that this is using the _mean_ of bootstrapped scores in results[[round]],
              # not the _median_ of bootstrapped scores in plots[[round]].
	      # This is for consistency with other plots (e.g., the heatmaps, where we
	      # report the means). We only use median in the boxplot of pearson (implicitly)
	      # and in the barplot of spearman (explicitly and for consistency with the pearson boxplot)
              mean.bootstrapped.scores <- results[[round]][["mean.bootstrapped.scores"]]
	      ## ... instead use median.bootstrapped.scores in plots[[round]]
	      ## median.bootstrapped.scores <- plots[[round]][["median.bootstrapped.scores"]]
              res <- llply(sub.challenges,   
                           .fun = function(sub.challenge) {
                               compute.chained.empirical.bayes(mean.bootstrapped.scores, bootstrapped.scores, sub.challenge, metric = "pearson")
                           })
          })


g.heatmap.merged.round1 <- plots[["1"]][["heatmaps"]][["merged"]]
## png(paste0(figs.dir, "fig-validation-round-1-merged-heatmap.png"), width = 2 * 480)
## print(g.heatmap.merged.round1)
## d <- dev.off()

g.strip.merged.round1.priority <- plots[["1"]][["strip.plots"]][["merged-priority"]]
g.strip.merged.round1.priority <- g.strip.merged.round1.priority +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 15),
          strip.text = element_text(size = 15))

## png(paste0(figs.dir, "fig-validation-round-1-merged-strip-plot-priority.png"), width = 2 * 480)
## print(g.strip.merged.round1.priority)
## d <- dev.off()

g <- plot_grid(g.strip.merged.round1.priority, g.heatmap.merged.round1, nrow = 2, labels = "AUTO")

png(paste0(figs.dir, "fig-validation-round-1-strip-and-heatmap-merged-cell-type", ".png"), width = 2 * 480, height = 2 * 480)                    
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-validation-round-1-strip-and-heatmap-merged-cell-type", ".pdf"), width = 2 * 7, height = 2 * 7)                    
print(g)
d <- dev.off()


cardinal.to.ordinal <- function(num) {
    lst <- list("1" = "First",
                "2" = "Up To Second",
                "3" = "Up To Third")
    if(!(num %in% names(lst))) {
        stop(paste0("Unknown ", num, "\n"))
    }
    return(lst[[num]])
}

##  - averaged heatmaps over rounds 2 and 3 (supp; one fig)
if(all(c("2", "3") %in% rounds)) {
  g.heatmap.merged.round2 <- plots[["2"]][["heatmaps"]][["merged"]]
  g.heatmap.merged.round2 <- g.heatmap.merged.round2 + ggtitle(paste0("Merged Coarse- and Fine-Grained (", cardinal.to.ordinal(2), " Submission)"))
  g.heatmap.merged.round2 <- g.heatmap.merged.round2 + theme(plot.title = element_text(hjust = 0.5))
  g.heatmap.merged.round3 <- plots[["3"]][["heatmaps"]][["merged"]]
  g.heatmap.merged.round3 <- g.heatmap.merged.round3 + ggtitle(paste0("Merged Coarse- and Fine-Grained (", cardinal.to.ordinal(3), " Submission)"))
  g.heatmap.merged.round3 <- g.heatmap.merged.round3 + theme(plot.title = element_text(hjust = 0.5))

  g <- plot_grid(g.heatmap.merged.round2, g.heatmap.merged.round3, nrow = 2, labels = "AUTO")

  png(paste0(figs.dir, "fig-validation-heatmap-rounds-2-and-3-merged-cell-type", ".png"), width = 2 * 480, height = 2 * 480)                    
  print(g)
  d <- dev.off()

  pdf(paste0(figs.dir, "fig-validation-heatmap-rounds-2-and-3-merged-cell-type", ".pdf"), width = 2 * 7, height = 2 * 7)                    
  print(g)
  d <- dev.off()
}

## - coarse and fine grained heatmap round 1
## - coarse and fine grained heatmap round 2
## - coarse and fine grained heatmap round 3
plot.rounds <- list("1" = "1", "2" = "2", "3" = "3")
plot.rounds <- plot.rounds[plot.rounds %in% rounds]
l_ply(plot.rounds,
      .fun = function(rd) {
          g.coarse <- plots[[rd]][["heatmaps"]][["coarse"]]
          g.coarse <- g.coarse + ggtitle(paste0("Coarse-Grained (", cardinal.to.ordinal(rd), " Submission)"))
          g.coarse <- g.coarse + theme(plot.title = element_text(hjust = 0.5))
          
          g.fine <- plots[[rd]][["heatmaps"]][["fine"]]
          g.fine <- g.fine + ggtitle(paste0("Fine-Grained (", cardinal.to.ordinal(rd), " Submission)"))
          g.fine <- g.fine + theme(plot.title = element_text(hjust = 0.5))
          
          g <- plot_grid(g.coarse, g.fine, nrow = 2, labels = "AUTO")
          png(paste0(figs.dir, "fig-validation-heatmap-round-", rd, "-coarse-and-fine-cell-type", ".png"), width = 2 * 480, height = 2 * 480)
          print(g)
          d <- dev.off()

          pdf(paste0(figs.dir, "fig-validation-heatmap-round-", rd, "-coarse-and-fine-cell-type", ".pdf"), width = 2 * 7, height = 2 * 7)
          print(g)
          d <- dev.off()
          
      })

##  - strip plots merged all methods rounds 1
##  - strip plots merged all methods rounds 2
##  - strip plots merged all methods rounds 3
plot.rounds <- list("1" = "1", "2" = "2", "3" = "3")
plot.rounds <- plot.rounds[plot.rounds %in% rounds]
l_ply(plot.rounds,
      .fun = function(rd) {
          g <- plots[[rd]][["strip.plots"]][["merged"]]
          g <- g + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 15),
                         strip.text = element_text(size = 15),
                         plot.title = element_text(hjust = 0.5))
          g <- g + ggtitle(paste0("Merged Coarse- and Fine-Grained (", cardinal.to.ordinal(rd), " Submission)"))
          
          png(paste0(figs.dir, "fig-validation-round-", rd, "-merged-strip-cell-type", ".png"), width = 2 * 480, height = 2 * 480)
          print(g)
          d <- dev.off()
          
          pdf(paste0(figs.dir, "fig-validation-round-", rd, "-merged-strip-cell-type", ".pdf"), width = 2 * 7, height = 2 * 7)
          print(g)
          d <- dev.off()

      })

##  - coarse and fine grained strip plots all methods round 1
##  - coarse and fine grained strip plots all methods round 2
##  - coarse and fine grained strip plots all methods round 3  
plot.rounds <- list("1" = "1", "2" = "2", "3" = "3")
plot.rounds <- plot.rounds[plot.rounds %in% rounds]
l_ply(plot.rounds,
      .fun = function(rd) {
          g1 <- plots[[rd]][["strip.plots"]][["coarse"]]
          g2 <- plots[[rd]][["strip.plots"]][["fine"]]          
          g1 <- g1 + theme(axis.text.y = element_text(size = 6),
                           axis.text.x = element_text(size = 15),
                           strip.text = element_text(size = 15),
                           plot.title = element_text(hjust = 0.5))
          g1 <- g1 + ggtitle(paste0("Coarse-Grained (", cardinal.to.ordinal(rd), " Submission)"))
          
          g2 <- g2 + theme(axis.text.y = element_text(size = 6),
                           axis.text.x = element_text(size = 15),
                           strip.text = element_text(size = 15),
                           plot.title = element_text(hjust = 0.5))
          g2 <- g2 + ggtitle(paste0("Fine-Grained (", cardinal.to.ordinal(rd), " Submission)"))         

          png(paste0(figs.dir, "fig-validation-round-", rd, "-coarse-and-fine-strip-cell-type", ".png"), width = 2 * 480, height = 2 * 480)
          g <- plot_grid(g1, g2, nrow = 2, labels = "AUTO")          
          print(g)
          d <- dev.off()

          pdf(paste0(figs.dir, "fig-validation-round-", rd, "-coarse-and-fine-strip-cell-type", ".pdf"), width = 2 * 7, height = 2 * 7)
          g <- plot_grid(g1, g2, nrow = 2, labels = "AUTO")          
          print(g)
          d <- dev.off()
          
      })

correlations <- list("spearman" = "spearman", "pearson" = "pearson")

apply.function <- function(results, fun = "min", ...) {
    res <-
        llply(results,
              .fun = function(res) {
                  llply(sub.challenges,
                        .fun = function(sub.challenge) {
                            llply(correlations, 
                                .fun = function(cor.method) {
                                    do.call(fun, list(res[["bootstrapped.scores"]][[sub.challenge]][, cor.method], ...))
                                }) }) }) }

x.mins <- apply.function(results, fun = "min", na.rm = TRUE)
x.maxs <- apply.function(results, fun = "max", na.rm = TRUE)

x.min.pearson.coarse <- x.mins[["1"]][["coarse"]][["pearson"]] 
x.min.pearson.coarse <- shift.limit(x.min.pearson.coarse)
x.max.pearson.coarse <- x.maxs[["1"]][["coarse"]][["pearson"]] 
x.max.pearson.coarse <- shift.limit(x.max.pearson.coarse)

x.min.pearson.fine <- x.mins[["1"]][["fine"]][["pearson"]] 
x.min.pearson.fine <- shift.limit(x.min.pearson.fine)
x.max.pearson.fine <- x.maxs[["1"]][["fine"]][["pearson"]] 
x.max.pearson.fine <- shift.limit(x.max.pearson.fine)
x.max.pearson.fine <- 0.85

x.min.spearman.coarse <- x.mins[["1"]][["coarse"]][["spearman"]] 
x.min.spearman.coarse <- shift.limit(x.min.spearman.coarse)
x.max.spearman.coarse <- x.maxs[["1"]][["coarse"]][["spearman"]] 
x.max.spearman.coarse <- shift.limit(x.max.spearman.coarse)

x.min.spearman.fine <- x.mins[["1"]][["fine"]][["spearman"]] 
x.min.spearman.fine <- shift.limit(x.min.spearman.fine)
x.max.spearman.fine <- x.maxs[["1"]][["fine"]][["spearman"]] 
x.max.spearman.fine <- shift.limit(x.max.spearman.fine)

g.bootstrap.coarse.pearson.round1 <- plots[["1"]][["barplots"]][["coarse-pearson"]] +
    scale_y_continuous(limits = c(x.min.pearson.coarse, x.max.pearson.coarse), expand = c(0, 0))
g.bootstrap.coarse.spearman.round1 <- plots[["1"]][["barplots"]][["coarse-spearman"]] +
    scale_y_continuous(limits = c(x.min.spearman.coarse, x.max.spearman.coarse), expand = c(0, 0))    

g.bootstrap.coarse.spearman.round1 <- g.bootstrap.coarse.spearman.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))
g.bootstrap.coarse.pearson.round1 <- g.bootstrap.coarse.pearson.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))

g.bootstrap.coarse.anno.round1 <- plots[["1"]][["barplots"]][["coarse-anno"]]
g.bootstrap.coarse.anno.legend.round1 <- plots[["1"]][["barplots"]][["coarse-legend"]]
g.bootstrap.coarse.anno.round1 <- g.bootstrap.coarse.anno.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))
g.bootstrap.coarse.anno.legend.round1 <- g.bootstrap.coarse.anno.legend.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))

g.bootstrap.coarse.pearson.round1 <- g.bootstrap.coarse.pearson.round1 + ylab("Pearson\nCorrelation") +
    theme(axis.title.y = element_blank())
g.bootstrap.coarse.spearman.round1 <- g.bootstrap.coarse.spearman.round1 + ylab("Spearman\nCorrelation")

title <- "Coarse-Grained (First Submission)"
## g.bootstrap.coarse.round1 <- grid.arrange(g.bootstrap.coarse.pearson.round1,
##                                          g.bootstrap.coarse.spearman.round1, nrow=1, widths = c(3,1),
##                                          top = textGrob(title, gp = gpar(fontsize = 20)))

plot_row <- plot_grid(g.bootstrap.coarse.pearson.round1,
                      g.bootstrap.coarse.spearman.round1, 
                      g.bootstrap.coarse.anno.round1,
                      g.bootstrap.coarse.anno.legend.round1, nrow=1, align="h", rel_widths = c(3,1.1,0.5,0.5))
g.bootstrap.coarse.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

g.bootstrap.fine.pearson.round1 <- plots[["1"]][["barplots"]][["fine-pearson"]] +
    scale_y_continuous(limits = c(x.min.pearson.fine, x.max.pearson.fine), expand = c(0, 0))
g.bootstrap.fine.spearman.round1 <- plots[["1"]][["barplots"]][["fine-spearman"]] +
    scale_y_continuous(limits = c(-0.05, x.max.spearman.fine), expand = c(0, 0))    

g.bootstrap.fine.spearman.round1 <- g.bootstrap.fine.spearman.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))
g.bootstrap.fine.pearson.round1 <- g.bootstrap.fine.pearson.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))

g.bootstrap.fine.pearson.round1 <- g.bootstrap.fine.pearson.round1 + ylab("Pearson\nCorrelation") +
    theme(axis.title.y = element_blank())
g.bootstrap.fine.spearman.round1 <- g.bootstrap.fine.spearman.round1 + ylab("Spearman\nCorrelation")

g.bootstrap.fine.anno.round1 <- plots[["1"]][["barplots"]][["fine-anno"]]
g.bootstrap.fine.anno.legend.round1 <- plots[["1"]][["barplots"]][["fine-legend"]]
g.bootstrap.fine.anno.round1 <- g.bootstrap.fine.anno.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))
g.bootstrap.fine.anno.legend.round1 <- g.bootstrap.fine.anno.legend.round1 +
    theme(text = element_text(size=18), title = element_text(size = 20))

title <- "Fine-Grained (First Submission)"
## g.bootstrap.fine.round1 <- grid.arrange(g.bootstrap.fine.pearson.round1,
##                                        g.bootstrap.fine.spearman.round1, nrow=1, widths = c(3,1),
##                                        top = textGrob(title, gp = gpar(fontsize = 20)))

plot_row <- plot_grid(g.bootstrap.fine.pearson.round1,
                      g.bootstrap.fine.spearman.round1, 
                      g.bootstrap.fine.anno.round1,
                      g.bootstrap.fine.anno.legend.round1, nrow=1, align="h", rel_widths = c(3,1.1,0.5,0.5))
g.bootstrap.fine.round1 <- plot_grid(textGrob(title, gp = gpar(fontsize = 20)), plot_row, ncol=1, rel_heights = c(0.1, 1))

g.round.coarse <- g.score.vs.round[["coarse"]]
g.round.coarse <- g.round.coarse + ggtitle("Coarse-Grained")
g.round.coarse <- g.round.coarse + theme(text = element_text(size=18), title = element_text(size = 20),
                                         axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=20))
g.round.coarse <- g.round.coarse + 
    scale_x_continuous(limits = c(x.min.pearson.coarse, x.max.pearson.coarse), expand = c(0, 0))
g.round.fine <- g.score.vs.round[["fine"]]
g.round.fine <- g.round.fine + ggtitle("Fine-Grained")
g.round.fine <- g.round.fine + theme(text = element_text(size=18), title = element_text(size = 20),
                                     axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=20))
g.round.fine <- g.round.fine +
    scale_x_continuous(limits = c(x.min.pearson.fine, x.max.pearson.fine), expand = c(0, 0))

g <- plot_grid(g.bootstrap.coarse.round1, g.bootstrap.fine.round1, g.round.coarse, g.round.fine, labels = c("A", "B", "C", "D"))

png(paste0(figs.dir, "fig-validation-round-1-performance", ".png"), width = 2 * 480, height = 2 * 480)
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-validation-round-1-performance", ".pdf"), width = 2 * 7, height = 2 * 7)
print(g)
d <- dev.off()



x.all.min.pearson.coarse <- min(unlist(lapply(x.mins, function(lst) lst[["coarse"]][["pearson"]])))
x.all.max.pearson.coarse <- min(unlist(lapply(x.maxs, function(lst) lst[["coarse"]][["pearson"]])))
x.all.min.spearman.coarse <- min(unlist(lapply(x.mins, function(lst) lst[["coarse"]][["spearman"]])))
x.all.max.spearman.coarse <- min(unlist(lapply(x.maxs, function(lst) lst[["coarse"]][["spearman"]])))

x.all.min.pearson.fine <- min(unlist(lapply(x.mins, function(lst) lst[["fine"]][["pearson"]])))
x.all.max.pearson.fine <- min(unlist(lapply(x.maxs, function(lst) lst[["fine"]][["pearson"]])))
x.all.min.spearman.fine <- min(unlist(lapply(x.mins, function(lst) lst[["fine"]][["spearman"]])))
x.all.max.spearman.fine <- min(unlist(lapply(x.maxs, function(lst) lst[["fine"]][["spearman"]])))

make.barplots <- function(plots, round,
                          x.min.pearson.coarse, x.max.pearson.coarse,
                          x.min.pearson.fine, x.max.pearson.fine,                          
                          x.min.spearman.coarse, x.max.spearman.coarse,
                          x.min.spearman.fine, x.max.spearman.fine) {                         

    g.bootstrap.coarse.pearson <- plots[[round]][["barplots"]][["coarse-pearson"]] +
        scale_y_continuous(limits = c(x.min.pearson.coarse, x.max.pearson.coarse), expand = c(0, 0))
    g.bootstrap.coarse.spearman <- plots[[round]][["barplots"]][["coarse-spearman"]] +
        scale_y_continuous(limits = c(x.min.spearman.coarse, x.max.spearman.coarse), expand = c(0, 0))    
    
    g.bootstrap.coarse.spearman <- g.bootstrap.coarse.spearman +
        theme(text = element_text(size=18), title = element_text(size = 20))
    g.bootstrap.coarse.pearson <- g.bootstrap.coarse.pearson +
        theme(text = element_text(size=18), title = element_text(size = 20))
    
    g.bootstrap.coarse.pearson <- g.bootstrap.coarse.pearson + ylab("Pearson\nCorrelation") +
        theme(axis.title.y = element_blank())
    g.bootstrap.coarse.spearman <- g.bootstrap.coarse.spearman + ylab("Spearman\nCorrelation")

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
        
    title <- paste0("Coarse-Grained (", sub.title, ")")
    g.bootstrap.coarse <- grid.arrange(g.bootstrap.coarse.pearson,
                                       g.bootstrap.coarse.spearman, nrow=1, widths = c(3,1),
                                       top = textGrob(title, gp = gpar(fontsize = 20)))
    
    g.bootstrap.fine.pearson <- plots[[round]][["barplots"]][["fine-pearson"]] +
        scale_y_continuous(limits = c(x.min.pearson.fine, x.max.pearson.fine), expand = c(0, 0))
    g.bootstrap.fine.spearman <- plots[[round]][["barplots"]][["fine-spearman"]] +
        scale_y_continuous(limits = c(x.min.spearman.fine, x.max.spearman.fine), expand = c(0, 0))    
    
    g.bootstrap.fine.spearman <- g.bootstrap.fine.spearman +
        theme(text = element_text(size=18), title = element_text(size = 20))
    g.bootstrap.fine.pearson <- g.bootstrap.fine.pearson +
        theme(text = element_text(size=18), title = element_text(size = 20))
    
    g.bootstrap.fine.pearson <- g.bootstrap.fine.pearson + ylab("Pearson\nCorrelation") +
        theme(axis.title.y = element_blank())
    g.bootstrap.fine.spearman <- g.bootstrap.fine.spearman + ylab("Spearman\nCorrelation")

    title <- paste0("Fine-Grained (", sub.title, ")")        
    g.bootstrap.fine <- grid.arrange(g.bootstrap.fine.pearson,
                                     g.bootstrap.fine.spearman, nrow=1, widths = c(3,1),
                                     top = textGrob(title, gp = gpar(fontsize = 20)))
    
    return(list("coarse" = g.bootstrap.coarse, "fine" = g.bootstrap.fine))
}


x.min.spearman.fine <- -0.05
g1.barplots <- make.barplots(plots, "1",
                             x.min.pearson.coarse, x.max.pearson.coarse, x.min.pearson.fine, x.max.pearson.fine,                          
                             x.min.spearman.coarse, x.max.spearman.coarse, x.min.spearman.fine, x.max.spearman.fine)

if("2" %in% plot.rounds) {
  g2.barplots <- make.barplots(plots, "2",
                               x.min.pearson.coarse, x.max.pearson.coarse, x.min.pearson.fine, x.max.pearson.fine,                          
                               x.min.spearman.coarse, x.max.spearman.coarse, x.min.spearman.fine, x.max.spearman.fine)
}

if("3" %in% plot.rounds) {
  g3.barplots <- make.barplots(plots, "3",
                               x.min.pearson.coarse, x.max.pearson.coarse, x.min.pearson.fine, x.max.pearson.fine,                          
                               x.min.spearman.coarse, x.max.spearman.coarse, x.min.spearman.fine, x.max.spearman.fine)
}

g <- plot_grid(
    g1.barplots$coarse, g1.barplots$fine,
    g2.barplots$coarse, g2.barplots$fine,
    g3.barplots$coarse, g3.barplots$fine, labels = "AUTO", nrow = 3)

png(paste0(figs.dir, "fig-validation-all-performance", ".png"), width = 2 * 480, height = 3 * 480)
print(g)
d <- dev.off()

pdf(paste0(figs.dir, "fig-validation-all-performance", ".pdf"), width = 2 * 7, height = 3 * 7)
print(g)
d <- dev.off()
