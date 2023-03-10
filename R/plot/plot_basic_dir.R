#setwd("/Users/daohanji/Desktop/Expression_Evolution/rerun/out_basic/dir/")
setwd("/Users/rexjiang/Desktop/Lab/multivariate_trait/expression/rerun/out_basic/dir")

# End-point result (directional selection)
dm<-read.table("1_1_out_basic_end.txt",sep="\t") # Data file: end-point phenotypes of replicate lineages that evolved under directional selection towards the same optimum

dnew=data.frame(dm[,1],dm[,2]) # Extract traits of interest
colnames(dnew)=c("mRNA","translation")
g<-ggplot(dnew,aes(x=mRNA,y=translation))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Translation rate") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor1-dir.pdf",plot=g,width=5.5,height=5) # Save the plot as a file

dnew=data.frame(dm[,1],dm[,3]) # Extract traits of interest
colnames(dnew)=c("mRNA","protein")
g<-ggplot(dnew,aes(x=mRNA,y=protein))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Protein level") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor2-dir.pdf",plot=g,width=5.5,height=5) # Save the plot as a file

# Use basic plot function
#plot(dm[,1],dm[,2],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
#abline(lm(dm[,2]~dm[,1]),col="blue",lwd=2)

#plot(dm[,1],dm[,3],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
#abline(lm(dm[,3]~dm[,1]),col="blue")

# End-point result (multiple optima)
dm<-read.table("multi_1_out_basic_end.txt",sep="\t") # Data file: end-point phenotypes of replicate lineages for all genes
opt.all=sort(unique(dm[,1])) # All optimal protein levels considered

# Colors to be used for different genes
color.all=c("darkgrey",
	"grey",
	"purple",
	"violet",
	"blueviolet",
	"blue",
	"lightblue",
	"aquamarine",
	"seagreen",
	"darkgreen",
	"green",
	"darkorange",
	"orange",
	"gold",
	"yellow",
	"red",
	"firebrick",
	"deeppink",
	"pink",
	"thistle")

colnames(dm)=c("opt","mRNA","translation","protein")
dm$opt=factor(dm$opt)

g<-ggplot(dm,aes(x=mRNA,y=translation,colour=opt))+geom_point(show.legend = FALSE)
g=g+scale_color_manual(values=color.all)
g=g+geom_smooth(method="lm",data=dm,linetype="dashed",colour="black") # Make a scatter plot and add the global regression line
# Add a regression line for each gene
for(i in 1:length(opt.all)){
	g=g+geom_smooth(method="lm",data=subset(dm,opt==opt.all[i]),colour=color.all[i])
}
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Translation rate") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor1-dir-multi.pdf",plot=g,width=5.5,height=5) # Save the plot as a file

g<-ggplot(dm,aes(x=mRNA,y=protein,colour=opt))+geom_point(show.legend = FALSE)
g=g+scale_color_manual(values=color.all)
g=g+geom_smooth(method="lm",data=dm,linetype="dashed",colour="black") # Make a scatter plot and add the global regression line
# Add a regression line for each gene
for(i in 1:length(opt.all)){
	g=g+geom_smooth(method="lm",data=subset(dm,opt==opt.all[i]),colour=color.all[i])
}
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Protein level") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor2-dir-multi.pdf",plot=g,width=5.5,height=5) # Save the plot as a file

# Use basic plot function
dsub=dm[which(dm[,1]==opt.all[1]),]
plot(dsub[,2],dsub[,3],xlim=c(-4.5,3),ylim=c(-4.5,3),xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,3])-0.1,max(dsub[,3])+0.1)
fit=lm(dsub[,3]~dsub[,2])
abline(fit,col=color.all[1],lwd=2)
for(i in 1:length(opt.all)){
	dsub=dm[which(dm[,1]==opt.all[i]),]
	par(new=TRUE)
	plot(dsub[,2],dsub[,3],xlim=c(-4.5,3),ylim=c(-4.5,3),xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
	clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,3])-0.1,max(dsub[,3])+0.1)
	fit=lm(dsub[,3]~dsub[,2])
	abline(fit,col=color.all[i],lwd=2)
}
clip(min(dm[,2])-0.5,max(dm[,2])+0.5,min(dm[,3])-0.5,max(dm[,3])+0.5)
abline(lm(dm[,3]~dm[,2]),col="grey",lwd=2,lty="dashed")

dsub=dm[which(dm[,1]==opt.all[1]),]
plot(dsub[,2],dsub[,4],xlim=c(-4.5,3),ylim=c(-6,3.5),xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,4])-0.1,max(dsub[,4])+0.1)
fit=lm(dsub[,4]~dsub[,2])
abline(fit,col=color.all[1],lwd=2)
for(i in 1:length(opt.all)){
	dsub=dm[which(dm[,1]==opt.all[i]),]
	par(new=TRUE)
	plot(dsub[,2],dsub[,4],xlim=c(-4.5,3),ylim=c(-6,3.5),xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
	clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,4])-0.1,max(dsub[,4])+0.1)
	fit=lm(dsub[,4]~dsub[,2])
	abline(fit,col=color.all[i],lwd=2)
}
clip(min(dm[,2])-0.5,max(dm[,2])+0.5,min(dm[,4])-0.5,max(dm[,4])+0.5)
abline(lm(dm[,4]~dm[,2]),col="grey",lwd=2,lty="dashed")


