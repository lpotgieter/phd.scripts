library(phytools)
library(phangorn)
setwd('/home/laura/Desktop/Laura_Ruppert/1_data_Laura/Cercospora/11_CAFE/
Cercospora_effectors/')
species_tree = read.newick(file = "Cercospora_maxlike_tree_root_midpoint.nwk",text)
force.ultrametric<-function(tree,method=c("nnls","extend")){
method<-method[1]
if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
rooted=TRUE,trace=0)
else if(method=="extend"){
h<-diag(vcv(tree))
d<-max(h)-h
ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
y=tree$edge[,2])
tree$edge.length[ii]<-tree$edge.length[ii]+d
} else
cat("method not recognized: returning input tree\n\n")
tree
}
species_tree_ultra = force.ultrametric(species_tree)
write.tree(species_tree_ultra, file = "tree_ultrametric.txt", append = FALSE, digits = 10,
tree.names = FALSE)
is.ultrametric(species_tree_ultra)
