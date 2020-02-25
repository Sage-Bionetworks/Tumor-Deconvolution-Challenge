use.synapse <- FALSE
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

## Patient-wise
plot.patient <- function(df, patient, x.val = "value.gt",
                           y.val = "value.cs",
                           xlab = "CyTOF (Normalized)",
                           ylab = "CIBERSORT (Normalized)",
                           show.identity = TRUE) {
    df <- subset(df, sample == patient)
    g <- ggplot(data = df, aes_string(x = x.val, y = y.val))
    g <- g + geom_point(aes(shape = variable, colour = sample))
    g <- g + scale_shape_manual("cell type", values=1:nlevels(df$variable))
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    g <- g + stat_smooth(method='lm')
    if(show.identity) {
        g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }
    cor <- as.numeric(round(cor.test(df[, x.val], df[, y.val])$estimate,digits=2))
    title <- paste0("Cell-type-wise: Pearson cor = ", cor)
    g <- g + ggtitle(title)
    g
}

g.cs.pt <- plot.patient(cs.gt.df, patient = "S1")

png("sdy364-cibersort-patient.png")
print(g.cs.pt)
d <- dev.off()

## Cell-type-wise
plot.cell.type <- function(df, cell.type, x.val = "value.gt",
                           y.val = "value.cs",
                           xlab = "CyTOF (Normalized)",
                           ylab = "CIBERSORT (Normalized)",
                           show.identity = TRUE,
                           show.rmse = FALSE,
                           show.correlation = TRUE) {
    df <- subset(df, variable == cell.type)
    g <- ggplot(data = df, aes_string(x = x.val, y = y.val))
    g <- g + geom_point(aes(shape = variable, colour = sample))
    g <- g + scale_shape_manual("cell type", values=1:nlevels(df$variable))
    g <- g + xlab(xlab)
    g <- g + ylab(ylab)
    if(show.correlation) {
        g <- g + stat_smooth(method='lm')
    }
    if(show.identity) {
        g <- g + geom_abline(slope = 1, intercept = 0, linetype = "dashed")
    }
    cor <- as.numeric(round(cor.test(df[, x.val], df[, y.val])$estimate,digits=2))
    title <- paste0("Cell-type-wise")
    if(show.correlation) { title <- paste0(title, ": Pearson cor = ", cor) }
    if(show.rmse) {
        rmse <- round(sqrt(sum((df[,x.val]-df[,y.val])^2)), digits=2)
        if(show.correlation) {
            title <- paste0(title, "; RMSE = ", rmse)
        } else {
            title <- paste0(title, ": RMSE = ", rmse)
        }
    }
    g <- g + ggtitle(title)
    g
}

g <- plot.cell.type(cs.gt.df, cell.type = "CD8_T_cells")

g <- plot.cell.type(cs.gt.df, cell.type = "Eosinophils")
png("sdy364-cibersort-eosinophils.png")
print(g)
d <- dev.off()


g.neutro.mcp <- plot.cell.type(mcp.gt.df, cell.type = "Neutrophils", y.val = "value.mcp",
                               ylab = "MCP-Counter", show.identity = FALSE)

g.neutro.cs <- plot.cell.type(cs.gt.df, cell.type = "Neutrophils")

library(gridExtra)
png("sdy364-neutrophils-cs-and-mcp.png", width = 2 * 480)
grid.arrange(g.neutro.cs, g.neutro.mcp, nrow = 2)
d <- dev.off()

png("sdy364-patient-and-cell-type-wise.png", width = 2 * 480)
grid.arrange(g.cs.pt, g.neutro.cs, nrow = 2)
d <- dev.off()

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
  
