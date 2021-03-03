library(LEA)
library(tidyverse)
### Example of analysis using snmf ###
# Creation of the genotype file: genotypes.geno.
# The data contain 400 SNPs for 50 individuals.
#data("tutorial")
#write.geno(tutorial.R, "genotypes.geno")
################
# running snmf #
################
project.snmf = snmf("yourfiles.GT.FORMAT.geno",
K = 1:20,
entropy = TRUE,
ploidy = 1,
repetitions = 10,
project = "new")
# plot cross-entropy criterion of all runs of the project
plot(project.snmf, cex = 1.2, col = "lightblue", pch = 19)
# get the cross-entropy of the 10 runs for K = 4
ce = cross.entropy(project.snmf, K = 4)
# select the run with the lowest cross-entropy for K = 4
best = which.min(ce)
#add names
metadata = read.table("your .d file", header =T)
indv_snmf = read_tsv("your isolate name filet", col_names = F)
names(indv_snmf) = "Sample"
datalist = list()
for (i in c(6, 7, 8, 9)){
best = which.min(cross.entropy(project.snmf, K = i))
temp = as.data.frame(Q(project.snmf, i, best))
temp= cbind(indv_snmf, temp)
temp = temp %>%
gather("Cluster", "Admix_coef", -"Sample") %>%
mutate(K=i)
datalist[[i]] = as.tibble(temp)
}
snmf_results_per_K = bind_rows(datalist)%>%
inner_join(., metadata, by = c("Sample" = "Isolate")) %>%
unite(ID, Country, Location, Host, col = "for_display", remove = F)
ggplot(snmf_results_per_K, aes(x = reorder(Sample, ID), y = Admix_coef, fill = Cluster,
text = for_display)) +
geom_bar(position = "stack", stat = "identity", show.legend = F) +
facet_grid(K~.) +
theme_bw() +
theme(axis.title = element_blank(),
axis.text.x = element_text(angle = 60, hjust = 1),
legend.title = element_blank())
# display the Q-matrix
#my.colors <- c("tomato", "lightblue",
"olivedrab", "gold")
#barchart(project.snmf, K = 4, run = best,
border = NA, space = 0, col = my.colors,
xlab = "Individuals", ylab = "Ancestry proportions",
main = "Ancestry matrix") -> bp
#axis(1, at = 1:length(bp$order),
labels = bp$order, las = 3, cex.axis = .4)
