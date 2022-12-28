setwd("/Users/daohanji/Desktop/Expression_Evolution")

library(phytools)

# Read the tree file
tr.text=read.table("tr_20.txt",sep="\t")
tr=read.tree(text=tr.text[1,1])

# Assign an arbitrary root edge such that the root can be clearly shown when the tree is plotted
tr$root.edge=max(tr$edge.length)/10

# Edges along which the optimal phenotype is different from the ancestral value are differently colored
color0=rep("black",nrow(tr$edge))
color1=rep("black",nrow(tr$edge));color1[38]="red"
color2=rep("black",nrow(tr$edge));color2[20:38]="red"

# Plot the tree with external branch(es) or clade(s) colored
colors=color0
plot(tr,edge.color=colors,show.tip.label=FALSE,edge.width=2,root.edge=TRUE)