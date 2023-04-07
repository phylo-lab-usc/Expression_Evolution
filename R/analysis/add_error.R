# Add errors to simulated phenotypes to evaluate error's effect
# Errors assumed to be normally distributed; distribution of absolute error size asumed to be the same for mRNA and protein levels

esd.all=c(.05,.04,.03,.02,.01) # SD of error
d<-read.table("out_basic_end_all.txt",sep="\t")
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

# Generate datasets with errors
setwd("/Users/daohanji/Desktop/Expression_Evolution/rerun/out_basic/error/")
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	dnew=matrix(0,nrow=nrow(d),ncol=ncol(d))
	dnew[,1]=d[,1]
	dnew[,2]=d[,2]
	dnew[,3]=d[,3]+rnorm(nrow(d),mean=0,sd=esd)
	dnew[,5]=d[,5]+rnorm(nrow(d),mean=0,sd=esd)
	dnew[,4]=dnew[,5]-dnew[,3]
	fn=paste("out_basic_end_error_",esd,".txt",sep="")
	write.table(dnew,file=fn,sep="\t")
}

# Calculate correlations and variances for datasets with errors
cor.all=list()
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	fn=paste("out_basic_end_error_",esd,".txt",sep="")
	d<-read.table(fn,sep="\t")
	colnames(d)=c("Ne","width","t1","t2","t3")
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
	fn=paste("out_basic_end_error_stat_",esd,".txt",sep="")
	write.table(dv,file=fn,sep="\t")
}

# Neutral
setwd("..")
d<-read.table("0_out_basic_end.txt",sep="\t")

setwd("/Users/daohanji/Desktop/Expression_Evolution/rerun/out_basic/error/")
cor.error.neutral=matrix(0,nrow=length(esd.all)+1,ncol=3)
cor.error.neutral[nrow(cor.error.neutral),2]=cor(d[,1],d[,2])
cor.error.neutral[nrow(cor.error.neutral),3]=cor(d[,1],d[,3])
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	dnew=matrix(0,nrow=nrow(d),ncol=ncol(d))
	dnew[,1]=d[,1]+rnorm(nrow(d),mean=0,sd=esd)
	dnew[,3]=d[,3]+rnorm(nrow(d),mean=0,sd=esd)
	dnew[,2]=dnew[,3]-dnew[,1]
	fn=paste("0_out_basic_end_error_",esd,".txt",sep="")
	write.table(dnew,file=fn,sep="\t")

	cor.error.neutral[n,1]=esd
	cor.error.neutral[n,2]=cor(dnew[,1],dnew[,2])
	cor.error.neutral[n,3]=cor(dnew[,1],dnew[,3])
}
write.table(cor.error.neutral,file="0_out_basic_end_error_stat.txt",sep="\t")

