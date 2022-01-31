## See https://stackoverflow.com/questions/27803710/ggplot2-divide-legend-into-two-columns-each-with-its-own-title
plot.anno.heatmap.with.multiple.legends <-
    function(df, id.col, anno.columns, anno.pals) {

        suppressPackageStartupMessages(p_load("RColorBrewer"))
        df <- df[, c(id.col, anno.columns)]

        ## Assume annotations are characters
        for(col in c(id.col, anno.columns)) {
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
