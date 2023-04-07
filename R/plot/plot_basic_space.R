# Generate heatmaps that show correlation between traits under different combinations of Ne and fitness function SD

library(ggplot2)

# Input file
d<-read.table("out_basic_end_all.txt",sep="\t") # Data file: variance through time under stabilizing selection
colnames(d)=c("Ne","width","t1","t2","t3") # Set column names

Ne.all=sort(unique(d[,1])) # All effective population sizes to consider
width.all=sort(unique(d[,2])) # All SDs of the Gaussian fitness function to consider
par.all=matrix(0,nrow=length(Ne.all)*length(width.all),ncol=2) # Data matrix containing all combinations of Ne and fitness function SD
row=0
for(i in 1:length(Ne.all)){
	for(j in 1:length(width.all)){
		row=row+1
		par.all[row,1]=Ne.all[i]
		par.all[row,2]=width.all[j]
	}
}

# Generate a data matrix that contains end-point correlations for each combination of Ne and fitness function SD
dv=matrix(0,nrow=nrow(par.all),ncol=7)
for(i in 1:nrow(par.all)){
	sub=which(d$Ne==par.all[i,1]&d$width==par.all[i,2])
	dsub=d[sub,]
	dv[i,1:2]=par.all[i,]
	dv[i,3]=var(dsub$t1)
	dv[i,4]=var(dsub$t2)
	dv[i,5]=var(dsub$t3)
	dv[i,6]=cor(dsub$t1,dsub$t2)
	dv[i,7]=cor(dsub$t1,dsub$t3)
}
dv=data.frame(dv)
colnames(dv)=c("Ne","width","var1","var2","var3","cor12","cor13")
write.table(dv,file="out_basic_end_stat.txt",sep="\t") # Save the data matrix

#dv<-read.table("out_basic_end_stat.txt",sep="\t")

# Make heatmaps
#dv$width=factor(dv$width,levels=c("100","50","10","5","1"),ordered=TRUE) # Convert the column to factor such that it is not treated as a continuous variable by ggplot2

dv$width=1/dv$width # Convert fitness function SD to its reciprocal such that strength of selection increases with the value
dv$width=factor(dv$width,levels=c("0.01","0.02","0.1","0.2","1"),ordered=TRUE)

h1<-ggplot(dv,aes(x=log10(Ne),y=width,fill=cor12))+geom_tile()
h1=h1+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL) # Set colors
h1=h1+theme_classic() # Clear background
h1=h1+xlab(expression(paste("Lo",g[10],N[e])))#+ylab("Strength of selection (reciprocal of fitness function SD)") # Use subscripts in axis labels
h1=h1+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend

h2<-ggplot(dv,aes(x=log10(Ne),y=width,fill=cor13))+geom_tile()
h2=h2+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL) # Set colors
h2=h2+theme_classic() # Clear background
h2=h2+xlab(expression(paste("Lo",g[10],N[e])))#+ylab("Strength of selection (reciprocal of fitness function SD)") # Use subscripts in axis labels
h2=h2+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust font size of axis labels, axis text, and legend

# Save the plots as files
ggsave("plot_basic_space_cor1.pdf",plot=h1,width=6,height=5)
ggsave("plot_basic_space_cor2.pdf",plot=h2,width=6,height=5)


