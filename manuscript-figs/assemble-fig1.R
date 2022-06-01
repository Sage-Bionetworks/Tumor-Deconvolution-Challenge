library(pacman)
p_load(png)
p_load(cowplot)
p_load(magick)
p_load(ggplot2)
#f1a <- readPNG("fig1a.png")
#f1b <- readPNG("fig1b.png")

# See https://stackoverflow.com/questions/46799022/how-to-add-an-in-memory-png-image-to-a-plot
#interpolate = TRUE
#m1a <- grid::rasterGrob(f1a, interpolate = interpolate)
#m1b <- grid::rasterGrob(f1b, interpolate = interpolate)

#g <- plot_grid(m1a,m1b,nrow=2,labels="AUTO")

p1 <- ggdraw() + draw_image("fig1a.png")
p2 <- ggdraw() + draw_image("fig1b.png")
g <- plot_grid(p1,p2,nrow=2,labels="AUTO")

ggsave("fig1.png", g, dpi=2400)

ggsave("fig1.pdf", g, dpi=2400)
