# Plot pairwise phenotypic divergence against pairwise divergence time

library(ggplot2)
library(phytools)

# Read the tree file
nsp=50 # Number of species
T=1e5
tfn=paste("tr_",nsp,".txt",sep="")
tr.text=read.table(tfn,sep="\t")
tr=read.tree(text=tr.text[1,1])
#tr$edge.length=tr$edge.length*T/vcv(tr)[1,1] # Rescale the tree such that tree height match the pre-set value

# Make a data matrix that contains all pairwise divergence times
dm=vcv(tr)[1,1]-vcv(tr) # Distance matrix (C matrix)
dist.all=matrix(0,nrow=(nsp^2-nsp)/2,ncol=3)
row=0
for(i in 1:nsp){
	for(j in 1:nsp){
		if(i<j){
			row=row+1
			dist.all[row,1]=i
			dist.all[row,2]=j
			dist.all[row,3]=dm[i,j]
		}
	}
}

# Process the dataset and convert tip phenotypes to pairwise divergence
width=1
fn=paste(width,"_out_tr",".txt",sep="")
d<-read.table(fn,sep="\t")
col.rna=c();col.translation=c();col.protein=c()
for(i in 1:nsp){
	col.rna=c(col.rna,3*i-2)
	col.translation=c(col.translation,3*i-1)
	col.protein=c(col.protein,3*i)
}
d.rna=d[,col.rna]
d.translation=d[,col.translation]
d.protein=d[,col.protein]

dnew.rna=matrix(0,nrow=nrow(dist.all),ncol=nrow(d))
dnew.translation=matrix(0,nrow=nrow(dist.all),ncol=nrow(d))
dnew.protein=matrix(0,nrow=nrow(dist.all),ncol=nrow(d))

for(i in 1:nrow(d)){
	for(np in 1:nrow(dist.all)){
		sp1=dist.all[np,1]
		sp2=dist.all[np,2]
		dnew.rna[np,i]=abs(d.rna[i,sp1]-d.rna[i,sp2])
		dnew.translation[np,i]=abs(d.translation[i,sp1]-d.translation[i,sp2])
		dnew.protein[np,i]=abs(d.protein[i,sp1]-d.protein[i,sp2])
	}
}

fn1=paste(width,"_out_tr_pw_mRNA.txt",sep="")
fn2=paste(width,"_out_tr_pw_translation.txt",sep="")
fn3=paste(width,"_out_tr_pw_protein.txt",sep="")
write.table(dnew.rna,file=fn1,sep="\t")
write.table(dnew.translation,file=fn2,sep="\t")
write.table(dnew.protein,file=fn3,sep="\t")

#dnew.rna<-read.table(fn1,sep="\t")
#dnew.translation<-read.table(fn2,sep="\t")
#dnew.protein<-read.table(fn3,sep="\t")

# Reorganize the data to prepare them for ggplot2
col1=rep(dist.all[,3],3) # Column of pairwise phylogenetic distance
col2=c(rep("mRNA",nrow(dist.all)),rep("Translation",nrow(dist.all)),rep("Protein",nrow(dist.all))) # Columns of traits
col3=c(rowMeans(dnew.rna),rowMeans(dnew.translation),rowMeans(dnew.protein)) # Columns of pairwise phenotypic divergence
dnew=data.frame(col1,col2,col3)
colnames(dnew)=c("t","Trait","div")
dnew$Trait=factor(dnew$Trait,levels=c("mRNA","Translation","Protein"),ordered=TRUE) # Order trait names to be shown in the color legend
g<-ggplot(dnew,aes(x=t,y=div,colour=Trait),pallete=c())+geom_point()+geom_smooth(method="loess") # Make scatter plot and add curves
g=g+theme_classic() # Clear background
g=g+xlab("Divergence time")+ylab("Phenotypic divergence") # Axis labels
g=g+scale_color_manual(values=c("orange","grey","purple")) # Color data points for different traits differently
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font size of axis labels, axis text, legend title, and legend text
ggsave("t-div-pw.pdf",plot=g,width=7.5,height=5) # Save the plot


