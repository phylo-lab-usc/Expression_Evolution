# Fit phenotypic data (given phylogeny) to BM or OU model

library(phytools)
library(geiger)

# Read the tree
nsp=50 # Number of species
tfn=paste("tr_",nsp,".txt",sep="")
tr.text=read.table(tfn,sep="\t")
tr=read.tree(text=tr.text[1,1])

# Read simulation result (end-point)
width=1
fn=paste(width,"_out_tr.txt",sep="")
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

# Calculate likelihood for each row
aic.rna=matrix(0,nrow=nrow(d),ncol=2)
aic.translation=matrix(0,nrow=nrow(d),ncol=2)
aic.protein=matrix(0,nrow=nrow(d),ncol=2)

for(i in 1:nrow(d)){
	x1=as.numeric(d.rna[i,]);names(x1)=as.character(1:nsp)
	x2=as.numeric(d.translation[i,]);names(x2)=as.character(1:nsp)
	x3=as.numeric(d.protein[i,]);names(x3)=as.character(1:nsp)
	
	f1a=fitContinuous(tr,x1,model="BM")
	f1b=fitContinuous(tr,x1,model="OU")
	f2a=fitContinuous(tr,x2,model="BM")
	f2b=fitContinuous(tr,x2,model="OU")
	f3a=fitContinuous(tr,x3,model="BM")
	f3b=fitContinuous(tr,x3,model="OU")
	
	aic.rna[i,1]=f1a$opt$aicc
	aic.translation[i,1]=f2a$opt$aicc
	aic.protein[i,1]=f3a$opt$aicc
	
	aic.rna[i,2]=f1b$opt$aicc
	aic.translation[i,2]=f2b$opt$aicc
	aic.protein[i,2]=f3b$opt$aicc
}
out=data.frame(aic.rna,aic.translation,aic.protein)
fn=paste(width,"_out_tr_aic.txt",sep="")
write.table(out,file=fn,sep="\t")

# Fraction of simulations that favor OU over BM
length(which(aic.rna[,2]<aic.rna[,1]))/nrow(d)
length(which(aic.translation[,2]<aic.translation[,1]))/nrow(d)
length(which(aic.protein[,2]<aic.protein[,1]))/nrow(d)

# Calculate AIC weights from AICs
width=0
fn=paste(width,"_out_tr_aic.txt",sep="")
d<-read.table(fn,sep="\t")
out=matrix(0,nrow=nrow(d),ncol=ncol(d))
for(i in 1:nrow(d)){
	aic.rna=as.numeric(d[i,1:2])
	w.rna=aicw(aic.rna)
	out[i,1]=w.rna[1,3]
	out[i,2]=w.rna[2,3]

	aic.translation=as.numeric(d[i,3:4])
	w.translation=aicw(aic.translation)
	out[i,3]=w.translation[1,3]
	out[i,4]=w.translation[2,3]

	aic.protein=as.numeric(d[i,5:6])
	w.protein=aicw(aic.protein)
	out[i,5]=w.protein[1,3]
	out[i,6]=w.protein[2,3]
}
fn=paste(width,"_out_tr_aic_weight.txt",sep="")
write.table(out,file=fn,sep="\t")

out.sum=matrix(0,nrow=length(esd.all)+1,ncol=7)
out.sum[nrow(out.sum),2:7]=colMeans(out)



# Analyze data subject to measurement error
esd.all=c(.05,.04,.03,.02,.01) # Error SDs to consider
for(n in 1:length(esd.all)){
	esd=esd.all[n]
	fn=paste(width,"_out_error_",esd,"_aic.txt",sep="")
	d<-read.table(fn,sep="\t")
	out=matrix(0,nrow=nrow(d),ncol=ncol(d))
	for(i in 1:nrow(d)){
		aic.rna=as.numeric(d[i,1:2])
		w.rna=aicw(aic.rna)
		out[i,1]=w.rna[1,3]
		out[i,2]=w.rna[2,3]

		aic.translation=as.numeric(d[i,3:4])
		w.translation=aicw(aic.translation)
		out[i,3]=w.translation[1,3]
		out[i,4]=w.translation[2,3]

		aic.protein=as.numeric(d[i,5:6])
		w.protein=aicw(aic.protein)
		out[i,5]=w.protein[1,3]
		out[i,6]=w.protein[2,3]
	}
	out.sum[n,1]=esd
	out.sum[n,2:7]=colMeans(out)
	fn=paste(width,"_out_error_",esd,"_aic_weight.txt",sep="")
	write.table(out,file=fn,sep="\t")
}

colnames(out.sum)=c("error_sd","wAICc_mRNA_BM","wAICc_mRNA_OU","wAICc_translation_BM","wAICc_translation_OU","wAICc_protein_BM","wAICc_protein_OU")
fn=paste(width,"_out_sum_error_tr_new.txt",sep="")
write.table(out.sum,file=fn,sep="\t")


