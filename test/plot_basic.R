setwd("/Users/daohanji/Desktop/Expression_Evolution/out_basic")

# Time-variance plot
d<-read.table("out_basic_var_all.txt",sep="\t")
dm=d[which(d[,1]==1e3&d[,2]==1),]

plot(dm[,3],dm[,4],xlab="T",ylab="Transcription rate (mRNA level) variance")
plot(dm[,3],dm[,5],xlab="T",ylab="Translation rate variance")
plot(dm[,3],dm[,6],xlab="T",ylab="Protein level variance")
plot(dm[1:200,3],dm[1:200,6],xlab="T",ylab="Protein level variance")
plot(dm[,3],dm[,6]/dm[,4],xlab="T",ylab="Variance ratio")

# Time-variance plot (neutral)
dm<-read.table("0_1e5_out_basic_var.txt",sep="\t")
plot(dm[,1],dm[,2],xlab="T",ylab="Transcription rate (mRNA level) variance")
plot(dm[,1],dm[,3],xlab="T",ylab="Translation rate variance")
plot(dm[,1],dm[,4],xlab="T",ylab="Protein level variance")
plot(dm[,1],dm[,4]/dm[,2],xlab="T",ylab="Variance ratio")

# End-point result
d<-read.table("out_basic_end_all.txt",sep="\t")
dm=d[which(d[,1]==1e3&d[,2]==1),]

plot(dm[,3],dm[,4],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
abline(lm(dm[,4]~dm[,3]),col="blue",lwd=2)

plot(dm[,3],dm[,5],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
abline(lm(dm[,5]~dm[,3]),col="blue")

# End-point result (neutral)
dm<-read.table("0_1e5_out_basic_end.txt",sep="\t")

plot(dm[,1],dm[,2],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
abline(lm(dm[,2]~dm[,1]),col="blue",lwd=2)

plot(dm[,1],dm[,3],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
abline(lm(dm[,3]~dm[,1]),col="blue")

# End-point result (directional selection)
dm<-read.table("1_1_1e+05_out_basic_end.txt",sep="\t")

plot(dm[,1],dm[,2],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
abline(lm(dm[,2]~dm[,1]),col="blue",lwd=2)

plot(dm[,1],dm[,3],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
abline(lm(dm[,3]~dm[,1]),col="blue")

# End-point result (multiple optima)
dm<-read.table("multi_1_1e+05_out_basic_end.txt",sep="\t")
plot(dm[,2],dm[,3],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
plot(dm[,2],dm[,4],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)

opt.all=sort(unique(dm[,1]))
color.all=c("darkgrey","purple","blue","lightblue","darkgreen","green","orange","yellow","red","pink")

dsub=dm[which(dm[,1]==opt.all[1]),]
plot(dsub[,2],dsub[,3],xlim=c(-4,2),ylim=c(-4,2),xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,3])-0.1,max(dsub[,3])+0.1)
fit=lm(dsub[,3]~dsub[,2])
abline(fit,col=color.all[1],lwd=2)
for(i in 1:length(opt.all)){
	dsub=dm[which(dm[,1]==opt.all[i]),]
	par(new=TRUE)
	plot(dsub[,2],dsub[,3],xlim=c(-4,2),ylim=c(-4,2),xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
	clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,3])-0.1,max(dsub[,3])+0.1)
	fit=lm(dsub[,3]~dsub[,2])
	abline(fit,col=color.all[i],lwd=2)
}
clip(min(dm[,2])-0.5,max(dm[,2])+0.5,min(dm[,3])-0.5,max(dm[,3])+0.5)
abline(lm(dm[,3]~dm[,2]),col="grey",lwd=2,lty="dashed")

dsub=dm[which(dm[,1]==opt.all[1]),]
plot(dsub[,2],dsub[,4],xlim=c(-5,2.5),ylim=c(-5,2.5),xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,4])-0.1,max(dsub[,4])+0.1)
fit=lm(dsub[,4]~dsub[,2])
abline(fit,col=color.all[1],lwd=2)
for(i in 1:length(opt.all)){
	dsub=dm[which(dm[,1]==opt.all[i]),]
	par(new=TRUE)
	plot(dsub[,2],dsub[,4],xlim=c(-5,2.5),ylim=c(-5,2.5),xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
	clip(min(dsub[,2])-0.1,max(dsub[,2])+0.1,min(dsub[,4])-0.1,max(dsub[,4])+0.1)
	fit=lm(dsub[,4]~dsub[,2])
	abline(fit,col=color.all[i],lwd=2)
}
clip(min(dm[,2])-0.5,max(dm[,2])+0.5,min(dm[,4])-0.5,max(dm[,4])+0.5)
abline(lm(dm[,4]~dm[,2]),col="grey",lwd=2,lty="dashed")


# Results for a broader parameter space
colnames(d)=c("Ne","width","t1","t2","t3")

Ne.all=sort(unique(d[,1])) # Effective population size
width.all=sort(unique(d[,2])) # SD of Gaussian fitness function
par.all=matrix(0,nrow=length(Ne.all)*length(width.all),ncol=2)
row=0
for(i in 1:length(Ne.all)){
	for(j in 1:length(width.all)){
		row=row+1
		par.all[row,1]=Ne.all[i]
		par.all[row,2]=width.all[j]
	}
}

dv=matrix(0,nrow=nrow(par.all),ncol=7)
for(i in 1:nrow(par.all)){
	sub=which(d$Ne==par.all[i,1]&d$width==par.all[i,2])
	dsub=d[sub,]
	dv[i,1:2]=par.all[i,]
	dv[i,3]=var(dsub$t1)
	dv[i,4]=var(dsub$t2)
	dv[i,5]=var(dsub$t3)
	dv[i,6]=cor(dsub$t1,dsub$t2)
	dv[i,7]=var(dsub$t3)/var(dsub$t1)
}
#dv[,1]=log10(dv[,1])
#dv[,2]=log10(dv[,2])
dv=data.frame(dv)
colnames(dv)=c("Ne","width","var1","var2","var3","cor12","ratio")
write.table(dv,file="out_basic_final.txt",sep="\t")
dv=dv[which(dv[,2]>=1),]

# Heatmap
dv$width=factor(dv$width,levels=c("1","5","10","50","100"),ordered=TRUE)

h1<-ggplot(dv,aes(x=log10(Ne),y=width,fill=cor12))+geom_tile()
h1=h1+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL)
h1=h1+theme_classic() # White background
h1=h1+xlab(expression(paste("Lo",g[10],N[e])))+ylab("Fitness function SD") # Use subscripts in axis labels
h1=h1+theme(axis.title=element_text(size=20),axis.text.x=element_text(size=15),,axis.text.y=element_text(size=15))

h2<-ggplot(dv,aes(x=log10(Ne),y=width,fill=ratio))+geom_tile()
h2=h2+scale_fill_gradient2(low="purple4",mid="white",high="darkorange",midpoint=1,limits=c(0,4),name=NULL)
h2=h2+theme_classic() # White background
h2=h2+xlab(expression(paste("Lo",g[10],N[e])))+ylab("Fitness function SD") # Use subscripts in axis labels
h2=h2+theme(axis.title=element_text(size=20),axis.text.x=element_text(size=15),,axis.text.y=element_text(size=15))

ggsave("plot_basic_space_cor.pdf",plot=h1,width=11,height=10)
ggsave("plot_basic_space_ratio.pdf",plot=h2,width=11,height=10)

