# Simulate coevolution of mRNA and protein levels along a tree
# Rescale the tree and examine how this could affect the observations
# Keep all simulation settings the same except the transformation

library(phytools)
library(geiger)

T=1e5 # Duration of each simulation/Tree height

# Fitness (univariate Guassion fitness function)
# Calculate fitness given distance to optimum (d) and SD of the Gaussian fitness function (a)
fitness <- function(d,a){
	if(a==0){ # Neutrality
		w=1
	}else{
		w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a)
	}
	return(w)
}

# Fixation probability
# Calculate fixation probability of a mutation given ancestral phenotype (x1), mutant phenotype (x2), SD of fitness function (a), and effective population size (Ne)
# x1, x2, and a are all numbers when a single trait is considered; they are vectors of the same length if multiple traits are considered
fix.prob <- function(x1,x2,a,Ne){
	wa=fitness(x1,a) # Ancestral fitness
	wm=fitness(x2,a) # Mutant fitness
	if(wa>0){
		s=(wm/wa)-1 # Coefficient of selection
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
		}
	}else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
		if(wm==0){ # Both ancestral and mutant fitness are close to 0
			p=1/(2*Ne) # The mutation is considered neutral
		}else{ # Ancestral fitness is close to 0 while mutant fitness isn't
			p=1 # The mutation is considered strongly beneficial
		}
	}
	return(p)
}

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

nsp=50 # Number of species
if.new.tr=0 # Whether to use an existing tree or generate a new tree
tfn=paste("tr_",nsp,".txt",sep="")
if(if.new.tr==1){
	tr=pbtree(b=1,d=0,n=nsp)
	tr$tip.label=1:nsp # Rename the tips such that tip and node labels are aligned
	height=vcv(tr)[1,1]
	tr$edge.length=tr$edge.length*T/height
	# Write the tree in parenthesis format for future reference
	write.tree(tr,file=tfn)
}else{
	tr.text=read.table(tfn,sep="\t")
	tr=read.tree(text=tr.text[1,1])
	tr$edge.length=tr$edge.length*T/vcv(tr)[1,1]
}

# Confirm the node labels
#plotTree(tr);tiplabels();nodelabels()

if.sp=rep(0,nrow(tr$edge))
for(i in 1:nrow(tr$edge)){
	if(length(which(tr$edge[,1]==tr$edge[i,2]))==0){
		if.sp[i]=1
	}
}
rsp=which(if.sp==1)

# Rescaling
par.all=c(0.5,0)
tr.all=list()
for(i in 1:length(par.all)){
	trrs=rescale(tr,model="lambda",lambda=par.all[i])
	tr.all[[i]]=trrs
}
#tr.original=tr
tr=NULL

Ne=1e3 # Effective population size
width=1 # SD of Gaussian fitness function

# Mutation parameters (assumed to be invariable among branches)
lambda.all=c(1,1) # Rates of mutations affecting transcription rate and (per-mRNA) translation rate
sig.all=c(0.1,0.1) # SD of effect size of mutations affecting transcription rate and (per-mRNA) translation rate

Ntest=500

col1=rep(0,length(par.all)*Ntest)
out.all=matrix(0,nrow=length(par.all)*Ntest,ncol=nsp*3)
for(n in 1:length(par.all)){
	tr=tr.all[[n]]
	for(test in 1:Ntest){
		v.end=matrix(0,nrow=nrow(tr$edge),ncol=3)
		for(nb in 1:nrow(tr$edge)){
			# Convert branch length to the number of time steps for the simulation; must be an integer
			Tb=as.integer(tr$edge.length[nb])
		
			# Matrix for genotypic values and phenotypes along the branch
			vb=matrix(0,nrow=3,ncol=(Tb+1))
		
			# Search for the ancestral branch to determine the starting genotypic value
			ances=which(tr$edge[,2]==tr$edge[nb,1])
			if(length(ances)==0){
				vb[,1]=c(0,0,0)
			}else{
				vb[,1]=v.end[ances,]
			}
		
			if(Tb>=1){
				# Simulate evolution (only if Tb is positive)
				for(t in 2:(Tb+1)){
					vb[,t]=vb[,(t-1)]
					nm=rpois(1,lambda=sum(lambda.all))
					if(nm>0){
						for(i in 1:nm){
							phe_ances=vb[3,t]
							type=sample(1:2,1,prob=lambda.all/sum(lambda.all)) # Determine whether the mutation affects transcription or translation
							effect=rnorm(1,mean=0,sd=sig.all[type]) # Determine mutation effect size
							phe_mutant=phe_ances+effect # Mutant protein level
							pf=fix.prob(phe_ances,phe_mutant,width,Ne) # Calculate fixation probability
							if.fix=rbinom(n=1,size=1,prob=pf) # Determine if the mutation gets fixed
							vb[type,t]=vb[type,t]+effect*if.fix # If the mutation if fixed, add its effect to the population mean phenotype
							vb[3,t]=vb[1,t]+vb[2,t]
						}
					}
				}
			}
			v.end[nb,]=vb[,(Tb+1)]
		}
		# Write endpoint phenotypes
		# Columns: species 1's mRNA level, species 1's translation rate, species 1's protein level, species 2's mRNA level...
		rw=(n-1)*Ntest+test
		col1[rw]=par.all[n]
		v.end.sp=v.end[rsp,]
		for(sp in 1:nsp){
			out.all[rw,(3*sp-2):(3*sp)]=v.end.sp[sp,]
		}
	}
}
# Output
out.all=cbind(col1,out.all)
write.table(out.all,file="out_tr_rescale.txt",sep="\t")

# Output file that contains correlations only
out.final=matrix(0,nrow=length(par.all),ncol=3)
for(n in 1:length(par.all)){
	dsub=out.all[which(out.all[,1]==par.all[n]),]
	trrs=tr.all[[n]]
	m.all=matrix(0,nrow=nrow(dsub),ncol=9)
	for(i in 1:nrow(dsub)){
		dss=as.numeric(dsub[i,2:ncol(dsub)])
		dnew=v2m(dss,nsp,3)
		rownames(dnew)=trrs$tip
		m=ratematrix(trrs,dnew)
		m.all[i,]=m2v(m)
	}
	ma=v2m(colMeans(m.all),3,3)
	mcor=cov2cor(ma)
	out.final[n,1]=par.all[n]
	out.final[n,2]=mcor[1,2]
	out.final[n,3]=mcor[1,3]
}
write.table(out.final,file="out_tr_rescale_final.txt")
#x=out[,1];y=out[,2]
#plot(x,y);lines(x[order(x)],y[order(y)],lwd=2)
