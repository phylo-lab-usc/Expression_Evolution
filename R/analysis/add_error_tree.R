# Add error to results of simulations along a tree
# Errors assumed to be normally distributed; distribution of absolute error size asumed to be the same for mRNA and protein levels

library(phytools)
library(geiger)

# Rearrange a matrix into a vector by concatenating all the rows
m2v <- function(m){
	v=rep(0,nrow(m)*ncol(m))
	k=0
	for(i in 1:nrow(m)){
		for(j in 1:ncol(m)){
			k=k+1
			v[k]=m[i,j]
		}
	}
	return(v)
}

# Rearrange a vector into a matrix given row and column numbers of the matrix
v2m <- function(v,r,c){
	m=matrix(0,nrow=r,ncol=c)
	k=1
	for(i in 1:r){
		m[i,]=v[k:(k+c-1)]
		k=k+c
	}
	return(m)
}

# Read the tree
nsp=50
tfn=paste("tr_",nsp,".txt",sep="")
tr.text=read.table(tfn,sep="\t")
tr=read.tree(text=tr.text[1,1])

# Read simulation result (end-point)
width=1
fn=paste(width,"_out_tr.txt",sep="")
d<-read.table(fn,sep="\t")

# Generate data with error
esd.all=c(.05,.04,.03,.02,.01)
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	dnew=matrix(0,nrow=nrow(d),ncol=ncol(d))
	for(i in 1:nsp){
		dnew[,3*i-2]=d[,3*i-2]+rnorm(nrow(d),mean=0,sd=esd)
		dnew[,3*i]=d[,3*i]+rnorm(nrow(d),mean=0,sd=esd)
		dnew[,3*i-1]=dnew[,3*i]-dnew[,3*i-2]
	}
	fn=paste(width,"_out_tr_error_",esd,".txt",sep="")
	write.table(dnew,file=fn,sep="\t")
}

# Calculate evolutionary correlation between traits for different error SDs
cor.error.all=matrix(0,nrow=length(esd.all)+1,ncol=2)
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	cor.error.all[n,1]=esd
	fn=paste(width,"_out_tr_error_",esd,".txt",sep="")
	d<-read.table(fn,sep="\t")
	mat.all=matrix(0,nrow=nrow(d),ncol=9)
	for(i in 1:nrow(d)){
		di=v2m(as.numeric(d[i,]),nsp,3)
		rownames(di)=1:nsp
		mat=ratematrix(tr,di)
		mat.all[i,]=m2v(mat)
	}
	mm=v2m(colMeans(mat.all),3,3)
	cor.error.all[n,2]=cov2cor(mm)[1,2]
}

# Read the R matrix calculated from error-free results
# The file is automatically read as an array, so conversion is needed
rmat.no.error<-read.table(paste(width,"_rmat_tr.txt",sep=""),sep="\t")
rmat.no.error=rbind(rmat.no.error[[1]],rmat.no.error[[2]],rmat.no.error[[3]])

# Data matrix containing correlations between traits for different error SDs
cor.error.all[nrow(cor.error.all),2]=cov2cor(rmat.no.error)[1,2]
write.table(cor.error.all,file="cor_error_out.txt",sep="\t")
