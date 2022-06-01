library(pacman)
p_load(png)
p_load(cowplot)
f1a <- readPNG("fig1a.png")
f1b <- readPNG("fig1b.png")

# See https://stackoverflow.com/questions/23807021/how-to-do-in-r-load-an-image-file-print-text-on-image-save-modified-image
m1a <- grid::rasterGrob(f1a, interpolate = TRUE)
m1b <- grid::rasterGrob(f1b, interpolate = TRUE)

g <- plot_grid(m1a,m1b,nrow=2,labels="AUTO")

png("fig1.png")
print(g)
d <- dev.off()

pdf("fig1.pdf")
print(g)
d <- dev.off()
