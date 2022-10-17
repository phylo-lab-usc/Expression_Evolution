# Modeling coevolution between mRNA and protein levels for two interacting genes

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
coeff=c(0.5,-0.5) # Interaction parameters; effect of gene 1 on gene 2 and effect of gene 2 on gene 1, respectively

# Fitness function
# Calculate fitness given distance to the optimal phenotype and shape of fitness function
fitness <- function(d,a){
	if(a==0){
		w=1
	}else{
		w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a) #Gaussian fitness function with SD equal to a
	}
	return(w)
}

# Fixation probability of a mutation
# Calculate fixation probability given ancestral phenotype, mutant phenotype, and shape of fitness function
fix.prob <- function(x1,x2,a){
	wa=fitness(x1,a) # Ancestral fitness
	wm=fitness(x2,a) # Mutant fitness
	if(wa>0){
		s=(wm/wa)-1 # Selection coefficient
		if(s==0){
			p=1/(2*Ne) # Neutral mutation
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s)) # Non-neutral mutation
		}
	}else{ # Ancestral fitness is recognized as 0 by R
		if(wm==0){ # If mutant fitness is also recognized as 0 by R, the mutation is considered neutral.
			p=1/(2*Ne)
		}else{ # If mutant fitness is not recognized as 0 by R, the mutation is considered strongly beneficial and would always fix.
			p=1
		}
	}
	return(p)
}

# Genotype-phenotype map
# Return the phenotype (equilibrium mRNA and protein abundances) given genotypic values and interaction parameters
# Input and output both in log scale
g2p <- function(gt,coeff){
	a1=gt[1];b1=gt[2];a2=gt[3];b2=gt[4]
	right=c(log(dcp)-b1,log(dcr)-a2,log(dcp)-b2,log(dcr)-a1)
	mat=rbind(c(1,-1,0,0),c(0,coeff[1],-1,0),c(0,0,1,-1),c(-1,0,0,coeff[2]))
	pt=solve(mat,right)
	return(pt) # Return log scale phenotypes
}

# Mutation rate of each trait (mean number of mutations per unit time)
# Equivalent to 2*Ne*u in population genetic models (u is mutation rate per genome per unit time)
lambda1=1 # Rate of mutations affecting transcription rate of gene 1
lambda2=1 # Rate of mutations affecting per-mRNA translation rate (translational efficiency) of gene 1
lambda3=1 # Rate of mutations affecting transcription rate of gene 2
lambda4=1 # Rate of mutations affecting per-mRNA translation rate (translational efficiency) of gene 2
lambda.all=lambda1+lambda2+lambda3+lambda4 # Total mutation rate

# Variance of mutation effect size (log scale)
sig1=0.1
sig2=0.1
sig3=0.1
sig4=0.1

width=1 # Width of fitness fucntion (x-axis in log scale)
T=1e4 # Duration of the simulation

# Matrix to store genotypic values
# Rows: mRNA level of gene 1, translational efficiency of gene 1, mRNA level of gene 2, translational efficiency of gene 2
gt=matrix(0,nrow=4,ncol=(T+1))
# Matrix to store phenotypes
# Rows: mRNA level of gene 1, protein level of gene 1, mRNA level of gene 2, protein level of gene 2
pt=matrix(0,nrow=4,ncol=(T+1))
pt[,1]=g2p(gt[,1],coeff) # Calculate initial phenotype

for(t in 2:(T+1)){
	gt[,t]=gt[,(t-1)]
	nm=rpois(1,lambda=lambda.all) # Total number of mutations to occur in this time step
	if(nm>0){
		for(i in 1:nm){
			gt_ances=gt[,t] # Ancestral genotype
			pt_ances=g2p(gt_ances,coeff) # Ancestral phenotype
			gt_mutant=gt_ances
			type=sample(1:4,1,prob=c(lambda1/lambda.all,lambda2/lambda.all,lambda3/lambda.all,lambda4/lambda.all)) # Decide which trait the mutation would affect
			if(type==1){
				effect=rnorm(1,mean=0,sd=sig1) # Mutation's phenotypic effect
				gt_mutant[1]=gt_mutant[1]+effect # Genotypic value of the mutant (original scale)
			}
			if(type==2){
				effect=rnorm(1,mean=0,sd=sig2) # Mutation's phenotypic effect
				gt_mutant[2]=gt_mutant[2]+effect # Genotypic value of the mutant (original scale)
			}
			if(type==3){
				effect=rnorm(1,mean=0,sd=sig3) # Mutation's phenotypic effect
				gt_mutant[3]=gt_mutant[3]+effect # Genotypic value of the mutant (original scale)
			}
			if(type==4){
				effect=rnorm(1,mean=0,sd=sig4) # Mutation's phenotypic effect
				gt_mutant[4]=gt_mutant[4]+effect # Genotypic value of the mutant (original scale)
			}
			pt_mutant=g2p(gt_mutant,coeff) # Mutant phenotype
			pf=fix.prob(pt_ances[4],pt_mutant[4],width) # Calculate fixation probability
			if.fix=rbinom(n=1,size=1,prob=pf) # Decide whether the mutation would fix
			if(if.fix==1){
				gt[,t]=gt_mutant # If the mutation fixes, the mutant genotype becomes the ancestral genotype when the next mutation would be examined
			}
		}
	}
	pt[,t]=g2p(gt[,t],coeff) # Record the phenotype (log scale)
}
#write.table(pt[,(T+1)],file="test_out.txt",sep="\t")

# Simulate multiple replicate lineages
Nrep=100 # Number of replicate lineages
gt=list() # Genotypic values of all lineages through time (each lineage would be a matrix in the list)
pt=list() # Phenotypes of all lineages through time (each lineage would be a matrix in the list)
for(n in 1:Nrep){
	gt[[n]]=matrix(0,nrow=4,ncol=(T+1)) # Genotypic values of the n-th lineage
	pt[[n]]=matrix(0,nrow=4,ncol=(T+1)) # Phenotypes of the n-th lineage
	pt[[n]][,1]=g2p(gt[[n]][,1],coeff) # Initial phenotype of the n-th lineage
	for(t in 2:(T+1)){
		gt[[n]][,t]=gt[[n]][,(t-1)]
		nm=rpois(1,lambda=lambda.all)
		if(nm>0){
			for(i in 1:nm){
				gt_ances=gt[[n]][,t]
				pt_ances=g2p(gt_ances,coeff)
				gt_mutant=gt_ances
				type=sample(1:4,1,prob=c(lambda1/lambda.all,lambda2/lambda.all,lambda3/lambda.all,lambda4/lambda.all))
				if(type==1){
					effect=rnorm(1,mean=0,sd=sig1)
					gt_mutant[1]=gt_mutant[1]+effect
				}
				if(type==2){
					effect=rnorm(1,mean=0,sd=sig2)
					gt_mutant[2]=gt_mutant[2]+effect
				}
				if(type==3){
					effect=rnorm(1,mean=0,sd=sig3)
					gt_mutant[3]=gt_mutant[3]+effect
				}
				if(type==4){
					effect=rnorm(1,mean=0,sd=sig4)
					gt_mutant[4]=gt_mutant[4]+effect
				}
				pt_mutant=g2p(gt_mutant,coeff)
				pf=fix.prob(pt_ances[4],pt_mutant[4],width)
				if.fix=rbinom(n=1,size=1,prob=pf)
				if(if.fix==1){
					gt[[n]][,t]=gt_mutant
				}
			}
		}
		pt[[n]][,t]=g2p(gt[[n]][,t],coeff)
	}
}
#write.table(pt[[Nrep]][,(T+1)],file="test_out.txt",sep="\t")

gt_end=matrix(0,nrow=Nrep,ncol=4) # Genotypic values of all lineages at the end of simulation
pt_end=matrix(0,nrow=Nrep,ncol=4) # Phenotypes of all lineages at the end of simulation
# Create a new directory to store output files
dir_new=paste(width,"_",coeff[1],"_",coeff[2])
dir.create(dir_new)
setwd(paste("./",dir_new,sep="")) # Change working directory to the newly created one
for(n in 1:Nrep){
	gt_end[n,]=gt[[n]][,(T+1)] # End-point genotypic valyes of this lineage
	pt_end[n,]=pt[[n]][,(T+1)] # End-point phenotype of the lineage

	# Create new directory for this lineage
	dir_new=paste("out_",n,sep="")
	dir.create(dir_new)
	fn1=paste(dir_new,"/gt_all.txt",sep="");write.table(gt[[n]],file=fn1,sep="\t") # Genotyic values through time for this lineage
	fn2=paste(dir_new,"/pt_all.txt",sep="");write.table(pt[[n]],file=fn2,sep="\t") # Phenotypes through time for this lineage
}
# Output end-point genotypic values and phenotypes
write.table(gt_end,file="gt_end.txt",sep="\t")
write.table(pt_end,file="pt_end.txt",sep="\t")
# Correlation matrix for all genotypic values and phenotypes
m.out=cov2cor(cov(data.frame(gt_end,pt_end)))
write.table(m.out,file="cor_mat.txt",sep="\t")

# Data matrix that contains variances of all traits through time; the last column contains end-point variances that would be used in most analyses
var.all=matrix(0,nrow=8,ncol=(T+1))
for(t in 2:(T+1)){
	d=matrix(0,nrow=Nrep,ncol=8)
	for(n in 1:Nrep){
		d[n,1:4]=gt[[n]][,t]
		d[n,5:8]=pt[[n]][,t]
	}
	for(i in 1:8){
		var.all[i,t]=var(d[,i])
	}
}
write.table(var.all,file="var_all.txt",sep="\t")

setwd("..")

