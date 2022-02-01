use.synapse <- TRUE
if(use.synapse) {
library(synapser)
synLogin()

sdy364.cibersort.results.synId <- "syn18345976"
sdy364.mcp.results.synId <- "syn18345978"
sdy364.gt.synId <- "syn18345202"

sdy364.cibersort.results.path <-
    synGet(sdy364.cibersort.results.synId)$path

sdy364.mcp.results.path <-
    synGet(sdy364.mcp.results.synId)$path

sdy364.gt.path <-
    synGet(sdy364.gt.synId)$path

} else {
    
    sdy364.cibersort.results.path <- "/Users/Brian/work/sage/deconvolution/Tumor-Deconvolution-Challenge/analysis/SDY364/cibersort_results_flow.tsv"
    
    sdy364.mcp.results.path <- "/Users/Brian/work/sage/deconvolution/Tumor-Deconvolution-Challenge/analysis/SDY364/mcpcounter_results_flow.tsv"
            
    sdy364.gt.path <- "/Users/Brian/.synapseCache/865/36847865/ground_truth_flow.tsv"
    
}
    

sdy364.cibersort.results.df <-
    read.table(sdy364.cibersort.results.path,
               sep="\t", header=TRUE)

sdy364.mcp.results.df <-
    read.table(sdy364.mcp.results.path,
               sep="\t", header=TRUE)
sdy364.mcp.results.df <- as.data.frame(t(sdy364.mcp.results.df))
sdy364.mcp.results.df$sample <- rownames(sdy364.mcp.results.df)
colnames(sdy364.mcp.results.df) <-
    gsub(colnames(sdy364.mcp.results.df), pattern=" ", replacement="_")

sdy364.gt.df <-
    read.table(sdy364.gt.path,
               sep="\t", header=TRUE)

apply.translation <- function(df, translation.lst) {
    for(i in 1:length(translation.lst)) {
        df[, names(translation.lst)[i]] <-
            unlist(apply(df[, translation.lst[[i]], drop = FALSE], 1,
                         function(row) { sum(as.numeric(row)) }))
    }
    df
}

gt.translation <- list(
    "Memory_B_cells" = c("IgD..memory.B.cell", "IgD..memory.B.cell.1"),
    "Eosinophils" = c("eosinophil"), 
    "CD4_T_cells" = c("CD4..T.cell"), 
    "CD8_T_cells" = c("CD8..T.cell"), 
    "Monocytes" = c("monocyte"), 
    "Myeloid_dendritic_cells" = c("myeloid.dendritic.cell"), 
    "Naive_B_cells" = c("na.ve.B.cell"), 
    "Neutrophils" = c("neutrophil"), 
    "T_cells" = c("T.cell"), 
    "T_follicular_helper_cell" = c("T.follicular.helper.cell"), 
    "B_cells" = c("B.cell")
)

sdy364.gt.df <-
    apply.translation(sdy364.gt.df, gt.translation)
sdy364.gt.df <-
    sdy364.gt.df[, c("sample", names(gt.translation))]
sdy364.gt.df$sample <- paste0("S", 1:nrow(sdy364.gt.df))

cs.translation <- list(
    "CD4_T_cells" = c("T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                      "T.cells.regulatory..Tregs.", "T.cells.gamma.delta"),
    "CD8_T_cells" = c("T.cells.CD8"),  
    "Eosinophils" = c("Eosinophils"),
    "Memory_B_cells" = c("B.cells.memory"),
    "Monocytes" = c("Monocytes"),
    "Myeloid_dendritic_cells" = c("Dendritic.cells.resting", "Dendritic.cells.activated"),
    "Naive_B_cells" = c("B.cells.naive"),
    "Neutrophils" = c("Neutrophils"),
    "T_cells" = c("T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                  "T.cells.regulatory..Tregs.", "T.cells.gamma.delta", "T.cells.CD8"),
    "T_follicular_helper_cells" = c("T.cells.regulatory..Tregs."),
    "B_cells" = c("B.cells.naive", "B.cells.memory")
)

sdy364.cibersort.results.df <-
    apply.translation(sdy364.cibersort.results.df, cs.translation)
sdy364.cibersort.results.df <-
    sdy364.cibersort.results.df[, c("sample", names(cs.translation))]
sdy364.cibersort.results.df$sample <- paste0("S", 1:nrow(sdy364.cibersort.results.df))

mcp.translation <- list(
    "CD8_T_cells" = c("CD8_T_cells"),
    "Monocytes" = c("Monocytic_lineage"),
    "Myeloid_dendritic_cells" = c("Myeloid_dendritic_cells"),
    "Neutrophils" = c("Neutrophils"),
    "T_cells" = c("T_cells"),
    "B_cells" = c("B_lineage")
)

sdy364.mcp.results.df <-
    apply.translation(sdy364.mcp.results.df, mcp.translation)
sdy364.mcp.results.df <-
    sdy364.mcp.results.df[, c("sample", names(mcp.translation))]
sdy364.mcp.results.df$sample <- paste0("S", 1:nrow(sdy364.mcp.results.df))

require(magrittr)
require(tibble)
require(plyr)
require(purrr)

df_to_matrix <- function(df, id_column){
    df %>% 
        data.frame() %>% 
        tibble::column_to_rownames(id_column) %>% 
        as.matrix()
}

matrix_to_df <- function(matrix, new_col){
    matrix %>% 
        data.frame() %>% 
        tibble::rownames_to_column(new_col) %>% 
        tibble::as_tibble()
}
library(reshape2)

cols <- intersect(colnames(sdy364.mcp.results.df),
                  colnames(sdy364.gt.df))
cols <- cols[!(cols == "sample")]

mcp.tmp <- sdy364.mcp.results.df[, c("sample", cols)]
gt.tmp <- sdy364.gt.df[, c("sample", cols)]

gt.tmp <-
    gt.tmp %>%
    df_to_matrix("sample") %>%
    divide_by(rowSums(.)) %>%
    matrix_to_df("sample")

mcp.melt <- melt(mcp.tmp, id = "sample")
gt.melt <- melt(gt.tmp, id = "sample")

mcp.gt.df <- merge(mcp.melt, gt.melt, by = c("sample", "variable"),
                   suffixes = c(".mcp", ".gt"))
mcp.gt.df$variable <- factor(mcp.gt.df$variable)

cols <- intersect(colnames(sdy364.cibersort.results.df),
                  colnames(sdy364.gt.df))
cols <- cols[!(cols == "sample")]

cs.tmp <- sdy364.cibersort.results.df[, c("sample", cols)]
gt.tmp <- sdy364.gt.df[, c("sample", cols)]

gt.tmp <-
    gt.tmp %>%
    df_to_matrix("sample") %>%
    divide_by(rowSums(.)) %>%
    matrix_to_df("sample")

cs.tmp <-
    cs.tmp %>%
    df_to_matrix("sample") %>%
    divide_by(rowSums(.)) %>%
    matrix_to_df("sample")

cs.melt <- melt(cs.tmp, id = "sample")
gt.melt <- melt(gt.tmp, id = "sample")

cs.gt.df <- merge(cs.melt, gt.melt, by = c("sample", "variable"),
                   suffixes = c(".cs", ".gt"))
cs.gt.df$variable <- factor(cs.gt.df$variable)

library(ggplot2)
plot.all <- function(df, x.val = "value.gt",
                     y.val = "value.cs",
                     xlab = "CyTOF (Normalized)",
                     ylab = "CIBERSORT (Normalized)",
                     show.identity = TRUE) {
    g <- ggplot(data = df, aes_string(x = x.val, y = y.val))
    g <- g + geom_point(aes(shape = variable, colour = sample))
    g <- g + scale_shape_manual("cell type", values=1:nlevels(df$variable))
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    g <- g + stat_smooth(method='lm')
    g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    cor <- as.numeric(round(cor.test(df[, x.val], df[, y.val])$estimate,digits=2))
    title <- paste0("Pearson cor = ", cor)
    g <- g + ggtitle(title)
    g
}

g.cs <- plot.all(cs.gt.df)
png("sdy364-cibersort-all.png")
print(g.cs)
d <- dev.off()

lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- ct$estimate
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    eq <- NULL
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=3, scientific=0),
                              pval = format(pval, digits=3, scientific=0)))
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=3, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}


add.correlation.text <- function(df, g, display.r2 = FALSE, method = "pearson", display.pval = FALSE, xoffset = 0.5, ...) {
  ylimits <- NULL
if(FALSE) {
  use.ggplot.2.2.1.limit.code <- TRUE
  if(use.ggplot.2.2.1.limit.code) {
    ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_ranges[[1]]$x.range
  } else {
    ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  }
}
  xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range

print(df)
  g <- g + geom_text(size = 10, x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]),
                     y =ylimits[1] + 0.8 * (ylimits[2] - ylimits[1]),
                     label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g
}


## Patient-wise
plot.patient <- function(df, patient, x.val = "value.gt",
                         y.val = "value.cs",
                         xlab = "CyTOF (Normalized)",
                         ylab = "CIBERSORT (Normalized)",
                         show.identity = TRUE, show.title = TRUE,
                         show.legend = TRUE) {
    df <- subset(df, sample == patient)
    g <- ggplot(data = df, aes_string(x = x.val, y = y.val))
    g <- g + geom_point(aes(shape = variable), colour = "black",
                        fill = "black", size = 5,
                        show.legend = show.legend)
    values <- 1:nlevels(df$variable)
    ## fill can only be used for shapes 21 to 25
    values <- seq(from=21,length.out=nlevels(df$variable))
    flag <- values > 25
    ## 25 is the max shape value; recycle to 1 if there are more
    values[flag] <- values[flag] - 25
    g <- g + scale_shape_manual("cell type", values=values)
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    g <- g + stat_smooth(method='lm')
    if(show.identity) {
        g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }
    cor <- as.numeric(round(cor.test(df[, x.val], df[, y.val])$estimate,digits=2))
    if(show.title) {
      title <- paste0("Cell-type-wise: Pearson cor = ", cor)
      g <- g + ggtitle(title)
    }
    new.df <- data.frame(x = df[,x.val], y = df[,y.val])
    g <- add.correlation.text(new.df, g)
    g <- g + theme(text = element_text(size=20))
    g
}

exclude <- c("Monocytes", "Eosinophils", "Myeloid_dendritic_cells",
             "Neutrophils")
g.cs.pt <- plot.patient(subset(cs.gt.df, !(variable %in% exclude)),
                        patient = "S1", show.title = FALSE,
                        show.legend = TRUE)

png("sdy364-cibersort-patient.png", width = 2 * 480)
print(g.cs.pt)
d <- dev.off()

## Cell-type-wise
plot.cell.type <- function(df, cell.type, x.val = "value.gt",
                           y.val = "value.cs",
                           xlab = "CyTOF (Normalized)",
                           ylab = "CIBERSORT (Normalized)",
                           show.identity = TRUE,
                           show.rmse = FALSE,
                           show.correlation = TRUE,
                           show.correlation.text = TRUE,
                           show.title = TRUE,
                           show.legend = TRUE) {
    df <- subset(df, variable == cell.type)
    g <- ggplot(data = df, aes_string(x = x.val, y = y.val))
    g <- g + geom_point(aes(shape = variable, colour = sample), size = 5,
                        show.legend = show.legend)
    ##    g <- g + scale_shape_manual("cell type", values=1:nlevels(df$variable))
    g <- g + scale_shape_manual("cell type", values=c(16))
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    if(show.correlation) {
        g <- g + stat_smooth(method='lm')
    }
    if(show.identity) {
        g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }
    cor <- as.numeric(round(cor.test(df[, x.val], df[, y.val], method="spearman")$estimate,digits=2))
    title <- paste0("Cell-type-wise")
    if(show.correlation) { title <- paste0(title, ": Spearman cor = ", cor) }
    if(show.rmse) {
        rmse <- round(sqrt(sum((df[,x.val]-df[,y.val])^2)), digits=2)
        if(show.correlation) {
            title <- paste0(title, "; RMSE = ", rmse)
        } else {
            title <- paste0(title, ": RMSE = ", rmse)
        }
    }
    if(show.title) {
        g <- g + ggtitle(title)
    }
    new.df <- data.frame(x = df[,x.val], y = df[,y.val])
    if(show.correlation.text) {
        g <- add.correlation.text(new.df, g, method = "pearson")
    }
    g <- g + theme(text = element_text(size=20))
    g
}

plot.cell.type.large  <- function(df, cell.type, x.val = "value.gt",
                           y.val = "value.cs",
                           xlab = "CyTOF (Normalized)",
                           ylab = "CIBERSORT (Normalized)",
                           show.rmse = FALSE) {
			   

  g <- plot.cell.type(df, cell.type = cell.type,
                      y.val = y.val,
                      ylab = ylab,
                      show.title = FALSE,
                      show.correlation.text = FALSE,
                      show.correlation = FALSE,
                      show.identity = FALSE,
                      show.legend = FALSE)
  g <- g + theme(text = element_text(size = 35),
                 axis.title.x = element_text(size = 35),
                 axis.title.y = element_text(size = 35))
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  g
}

g <- plot.cell.type(cs.gt.df, cell.type = "CD8_T_cells")

g <- plot.cell.type(cs.gt.df, cell.type = "Eosinophils")
png("sdy364-cibersort-eosinophils.png")
print(g)
d <- dev.off()



g.b.cs <- plot.cell.type.large(cs.gt.df, cell.type = "B_cells",
                                y.val = "value.cs",
                                ylab = "CIBESORT Score")
g.b.cs <- g.b.cs + xlab("Ground Truth") + ylab("Prediction")
png("sdy364-bcell-cs-large.png")
print(g.b.cs)
d <- dev.off()

g.mono.cs <- plot.cell.type.large(cs.gt.df, cell.type = "Monocytes",
                                y.val = "value.cs",
                                ylab = "CIBESORT Score")
g.mono.cs <- g.mono.cs + xlab("Ground Truth") + ylab("Prediction")
png("sdy364-mono-cs-large.png")
print(g.mono.cs)
d <- dev.off()


g.neutro.mcp <- plot.cell.type(mcp.gt.df, cell.type = "Neutrophils", y.val = "value.mcp",
                               ylab = "MCP-Counter", show.identity = FALSE)

g.neutro.cs <- plot.cell.type(cs.gt.df, cell.type = "Neutrophils", show.title = FALSE)

g.cd8.cs <- plot.cell.type(cs.gt.df, cell.type = "CD8_T_cells",
                           show.title = FALSE,
                           show.correlation.text = FALSE,
                           show.legend = TRUE)

ct <- "CD8_T_cells"
ct <- "B_cells"
g.b.mcp <- plot.cell.type(mcp.gt.df, cell.type = ct,
                            y.val = "value.mcp",
                            ylab = "MCP-Counter Score",
                            show.title = FALSE,
                            show.correlation.text = FALSE,
                            show.correlation = FALSE,
                            show.identity = FALSE,
                            show.legend = TRUE)


png("sdy364-cd8-cs.png", width = 2 * 480)
print(g.cd8.cs)
d <- dev.off()

png("sdy364-bcell-mcp.png", width = 2 * 480)
print(g.b.mcp)
d <- dev.off()

if(FALSE) {
g.b.mcp <- plot.cell.type(mcp.gt.df, cell.type = ct,
                            y.val = "value.mcp",
                            ylab = "MCP-Counter Score",
                            show.title = FALSE,
                            show.correlation.text = FALSE,
                            show.correlation = FALSE,
                            show.identity = FALSE,
                            show.legend = FALSE)
g.b.mcp <- g.b.mcp + theme(text = element_text(size = 35),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35))
g.b.mcp <- g.b.mcp + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
g.b.mcp <- plot.cell.type.large(mcp.gt.df, cell.type = ct,
                                y.val = "value.mcp",
                                ylab = "MCP-Counter Score")
g.b.mcp <- g.b.mcp + xlab("Ground Truth") + ylab("Prediction")
png("sdy364-bcell-mcp-large.png")
print(g.b.mcp)
d <- dev.off()

g.mono.mcp <- plot.cell.type.large(mcp.gt.df, cell.type = "Monocytes",
                                y.val = "value.mcp",
                                ylab = "MCP-Counter Score")
g.mono.mcp <- g.mono.mcp + xlab("Ground Truth") + ylab("Prediction")
png("sdy364-mono-mcp-large.png")
print(g.mono.mcp)
d <- dev.off()



g.b.mcp <- plot.cell.type(mcp.gt.df, cell.type = ct,
                            y.val = "value.mcp",
                            ylab = "MCP-Counter Score",
                            show.title = FALSE,
                            show.correlation.text = TRUE,
                            show.correlation = TRUE,
                            show.identity = FALSE,
                            show.legend = TRUE)


png("sdy364-cd8-cs-corr.png", width = 2 * 480)
print(g.cd8.cs)
d <- dev.off()

png("sdy364-bcell-mcp-corr.png", width = 2 * 480)
print(g.b.mcp)
d <- dev.off()


library(gridExtra)
png("sdy364-neutrophils-cs-and-mcp.png", width = 2 * 480)
grid.arrange(g.neutro.cs, g.neutro.mcp, nrow = 2)
d <- dev.off()

library(patchwork)
##png("sdy364-patient-and-cell-type-wise.png", width = 2 * 480)
g.cs.pt <- g.cs.pt + labs(subtitle = "A")
g.cd8.cs <- g.cd8.cs + labs(subtitle = "B")
g.cs.pt + g.cd8.cs + plot_layout(nrow=2)
ggsave("sdy364-patient-and-cell-type-wise.png")
## grid.arrange(g.cs.pt, g.neutro.cs, nrow = 2)
##d <- dev.off()



ct <- "Myeloid_dendritic_cells"
g.mcp <- plot.cell.type(mcp.gt.df, cell.type = ct, y.val = "value.mcp",
                        ylab = "MCP-Counter", show.identity = FALSE,
                        show.rmse = TRUE)
g.cs <- plot.cell.type(cs.gt.df, cell.type = ct, show.rmse = TRUE)

png("sdy364-dc-cs-and-mcp.png", width = 2 * 480)
grid.arrange(g.cs, g.mcp, nrow = 2)
d <- dev.off()

mcp.gt.rescaled.df <- mcp.gt.df
mcp.gt.rescaled.df$value.gt <- mcp.gt.rescaled.df$value.gt /
    max(mcp.gt.rescaled.df$value.gt)
mcp.gt.rescaled.df$value.mcp <- mcp.gt.rescaled.df$value.mcp /
    max(mcp.gt.rescaled.df$value.mcp)

cs.gt.rescaled.df <- cs.gt.df
cs.gt.rescaled.df$value.gt <- cs.gt.rescaled.df$value.gt /
    max(cs.gt.rescaled.df$value.gt)
cs.gt.rescaled.df$value.cs <- cs.gt.rescaled.df$value.cs /
    max(cs.gt.rescaled.df$value.cs)

ct <- "Myeloid_dendritic_cells"
g.mcp <- plot.cell.type(mcp.gt.rescaled.df, cell.type = ct, y.val = "value.mcp",
                        ylab = "MCP-Counter", show.identity = FALSE, show.rmse = TRUE,
                        show.correlation = TRUE)
g.cs <- plot.cell.type(cs.gt.rescaled.df, cell.type = ct, show.rmse = TRUE,
                       show.identity = FALSE, show.correlation = TRUE)


png("sdy364-dc-cs-and-mcp-rescaled.png", width = 2 * 480)
grid.arrange(g.cs, g.mcp, nrow = 2)
d <- dev.off()

##- sdy364
##  - slide describing transformation to [0:1] across all data
##  - plots showing score vs transformed score
  
