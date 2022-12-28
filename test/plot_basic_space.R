setwd("/Users/daohanji/Desktop/Expression_Evolution")

library(ggplot2)

# Input file
d<-read.table("out_basic_end_all.txt",sep="\t")
# Set column names
colnames(d)=c("Ne","width","t1","t2","t3")

Ne.all=c(1e2,1e3,1e4,1e5) # Effective population size
width.all=c(0.1,1,10,100) # SD of Gaussian fitness function
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
dv[,1]=log10(dv[,1])
dv[,2]=log10(dv[,2])
dv=data.frame(dv)
colnames(dv)=c("Ne","width","var1","var2","var3","cor12","ratio")

h1<-ggplot(dv,aes(x=Ne,y=width,fill=cor12))+geom_tile()
h1=h1+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL)
h1=h1+theme_classic() # White background
h1=h1+xlab(bquote(N[e]))+ylab("Fitness function SD") # Use subscripts in axis labels

h2<-ggplot(dv,aes(x=Ne,y=width,fill=ratio))+geom_tile()
h2=h2+scale_fill_gradient2(low="purple4",mid="white",high="darkorange",midpoint=1,limits=c(0,4),name=NULL)
h2=h2+theme_classic() # White background
h2=h2+xlab(bquote(N[e]))+ylab("Fitness function SD") # Use subscripts in axis labels

ggsave("plot_basic_space_cor.pdf",plot=h1,width=11,height=10)
ggsave("plot_basic_space_ratio.pdf",plot=h2,width=11,height=10)
