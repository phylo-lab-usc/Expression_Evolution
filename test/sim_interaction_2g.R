# Modeling coevolution between mRNA and protein levels for two interacting genes
# Only one gene's protein level is directly under selection

setwd("/Users/daohanji/Desktop/Expression_Evolution")

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
coeff=c(0.5,-0.5) # Interaction parameters; effect of gene 1 on gene 2 and effect of gene 2 on gene 1, respectively

# Fitness function
# Calculate fitness given distance to optimum and SD of Gaussian fitness function
fitness <- function(d,a){ 
	if(a==0){
		w=1 # Fitness is set to be constant when there is no selection (a=0)
	}else{
		w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a) # Calculate fitness (normalized, between 0 and 1)
	}
	return(w)
}

# Fixation probability
# Calculate fixation probability given ancestral and mutant phenotypes (distances to optimum) and SD of Gaussian fitness function
fix.prob <- function(x1,x2,a){
	wa=fitness(x1,a) # Calculate ancestral fitness from ancestral phenotype
	wm=fitness(x2,a) # Calculate mutant fitness from mutant phenotype
	if(wa>0){
		s=(wm/wa)-1 # Coefficient of selection
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s)) # Fixation probability
		}
	}else{
		if(wm==0){ # If mutant fitness is also recognized as 0 by R, the mutation is considered neutral.
			p=1/(2*Ne)
		}else{ # If mutant fitness is not recognized as 0 by R, the mutation is considered strongly beneficial and would always fix.
			p=1
		}
	}
	return(p)
}

# Genotype-phenotype map
# Calculate the phenotype (equilibrium mRNA and protein abundances) given genotypic values and interaction parameters
# Input and output are both in log scale
g2p <- function(gt,coeff){
	a1=gt[1];b1=gt[2];a2=gt[3];b2=gt[4]
	right=c(log(dcp)-b1,log(dcr)-a2,log(dcp)-b2,log(dcr)-a1)
	mat=rbind(c(1,-1,0,0),c(0,coeff[1],-1,0),c(0,0,1,-1),c(-1,0,0,coeff[2]))
	pt=solve(mat,right)
	return(pt) # Return log scale phenotypes
}

# Mutation rate of each trait (mean number of mutations per unit time)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
lambda.all=c(1,1,1,1)

# Variance of mutation effect size (log scale)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
sig.all=c(0.1,0.1,0.1,0.1)

width=1 # Width of fitness fucntion (x-axis in log scale)
T=1e4 # Duration of the simulation

# Matrix to store genotypic values through time
# Rows: mRNA level of gene 1, translational efficiency of gene 1, mRNA level of gene 2, translational efficiency of gene 2
gt=matrix(0,nrow=4,ncol=(T+1))

# Matrix to store phenotypes through time
# Rows: mRNA level of gene 1, protein level of gene 1, mRNA level of gene 2, protein level of gene 2
pt=matrix(0,nrow=4,ncol=(T+1))
pt[,1]=g2p(gt[,1],coeff) # Calculate initial phenotype

# Simulate multiple replicate lineages
Nrep=100 # Number of replicate lineages
gt=list() # Genotypic values of all lineages through time (each lineage would be a matrix in the list)
pt=list() # Phenotypes of all lineages through time (each lineage would be a matrix in the list)
gt_end=matrix(0,nrow=Nrep,ncol=4) # Genotypic values of all lineages at the end of simulation
pt_end=matrix(0,nrow=Nrep,ncol=4) # Phenotypes of all lineages at the end of simulation
for(n in 1:Nrep){
	gt[[n]]=matrix(0,nrow=4,ncol=(T+1)) # Genotypic values of the n-th lineage
	pt[[n]]=matrix(0,nrow=4,ncol=(T+1)) # Phenotypes of the n-th lineage
	pt[[n]][,1]=g2p(gt[[n]][,1],coeff) # Initial phenotype of the n-th lineage
	for(t in 2:(T+1)){
		gt[[n]][,t]=gt[[n]][,(t-1)]
		nm=rpois(1,lambda=sum(lambda.all))
		if(nm>0){
			for(i in 1:nm){
				gt_ances=gt[[n]][,t] # Ancestral genotypic values
				pt_ances=g2p(gt_ances,coeff) # Ancestral phenotype
				gt_mutant=gt_ances
				type=sample(1:4,1,prob=lambda.all/sum(lambda.all))
				effect=rnorm(1,mean=0,sd=sig.all[type])
				gt_mutant[type]=gt_mutant[type]+effect # Mutant genotypic values
				pt_mutant=g2p(gt_mutant,coeff) # Mutant phenotype
				pf=fix.prob(pt_ances[4],pt_mutant[4],width)
				if.fix=rbinom(n=1,size=1,prob=pf)
				if(if.fix==1){
					gt[[n]][,t]=gt_mutant # The mutant genotypic values becomes the ancestral genotypic values before the next mutation would be considered
				}
			}
		}
		pt[[n]][,t]=g2p(gt[[n]][,t],coeff) # Phenotype of lineage n at time t
	}
	gt_end[n,]=gt[[n]][,(T+1)] # End-point genotypic valyes of this lineage
	pt_end[n,]=pt[[n]][,(T+1)] # End-point phenotype of the lineage
}

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

# Create a new directory to store output files
dir_new=paste(width,"_",coeff[1],"_",coeff[2])
dir.create(dir_new)
setwd(paste("./",dir_new,sep="")) # Change working directory to the newly created one
write.table(gt_end,file="gt_end.txt",sep="\t")
write.table(pt_end,file="pt_end.txt",sep="\t")
# Correlation matrix for all genotypic values and phenotypes
m.out=cov2cor(cov(data.frame(gt_end,pt_end)))
write.table(m.out,file="cor_mat.txt",sep="\t")
write.table(var.all,file="var_all.txt",sep="\t")
setwd("..")
