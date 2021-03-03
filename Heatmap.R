#heatmap in R:
#plot_cazymes.R
library(ggplot2)
data_h <- read.table("cazy_all_species_heatmap.txt", header = T)
mine.heatmap <- ggplot(data_h, mapping = aes(x = species, y = cazyme_family, fill = frequency))
+
geom_tile() +
theme(axis.text.y = element_text(size=6)) +
xlab(label = "Cercospora species") +
ylab(label = "CAZyme family") +
scale_fill_gradient(low = "grey87",
high = "grey20")
mine.heatmap
