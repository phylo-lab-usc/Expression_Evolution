# Modeling coevolution between mRNA and protein levels for three interacting genes

library(gtools)

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
ngene=3 # Number of genes under consideration

# Interaction parameters
# All permutations of 0, a positive number and a negative number
# Numbers represent elements [1,2], [1,3], [2,1], [2,3], [3,1], [3,2], respectively
permut.all=permutations(3,6,c(-0.5,0,0.5),repeats.allowed=TRUE)

# Get transpose of 3x3 matrix given off-diagonal elements
# "tpv" for "'transpose' of vector"
tpv <- function(v){
	new=rep(0,6)
	new[1]=v[3];new[3]=v[1]
	new[2]=v[5];new[5]=v[2]
	new[4]=v[6];new[6]=v[4]
	return(new)
}
# Remove redundant combinations
coeff.all=list()
num=0
for(i in 1:nrow(permut.all)){
	check=list(tpv(permut.all[i,])) %in% coeff.all
	if(check==FALSE){
		num=num+1
		coeff.all[[num]]=permut.all[i,]
	}
}

# Convert to matrix format
coeff.all.m=list()
for(i in 1:length(coeff.all)){
	m=matrix(0,nrow=3,ncol=3)
	m[1,2:3]=coeff.all[[i]][1:2]
	m[2,1]=coeff.all[[i]][3];m[2,3]=coeff.all[[i]][4]
	m[3,1:2]=coeff.all[[i]][5:6]
	coeff.all.m[[i]]=m
}

# Fitness function
# Calculate fitness given distance to the optimal phenotype and shape of fitness function
fitness <- function(d,a){
	w=rep(0,length(d))
	for(i in 1:length(d)){
		if(a[i]!=0){
			w[i]=dnorm(d[i],mean=0,sd=a[i])/dnorm(0,mean=0,sd=a[i]) # Calculate a fitness value for each trait
		}
	}
	w=w[which(a!=0)] # Remove traits that do not affect fitness
	w=prod(w)^(1/length(w)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
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
	ngene=length(gt)/2
	right=rep(0,length(gt))
	mat=matrix(0,nrow=length(gt),ncol=length(gt))
	for(i in 1:ngene){
		right[2*i-1]=log(dcr)-gt[2*i-1] # Element 2*i-1 of the right-hand term, B
		right[2*i]=log(dcp)-gt[2*i] # Element 2*i of the right-hand term, B
		for(j in 1:ngene){
			# If the row number is odd and the column number is even, the entry is an interaction parameter
			mat[2*i-1,2*j]=coeff[j,i]
			# If both the row number and the column number are odd, the entry is zero
		}
		mat[2*i-1,2*i-1]=-1 # Diagonal elements are all -1
		mat[2*i,2*i-1]=1 # If the row number is even and the column number is odd, the entry is 1
		mat[2*i,2*i]=-1 # Diagonal elements are all -1
	}
	pt=solve(mat,right) # Solve Ax=B equation system
	return(pt) # Return log scale phenotypes
}

# Rows of the output data that correspond to a certain type of phenotype
# Used to conveniently extract subset of data during simulations and analyses
row.transcription=c() # Rows of output data that are transcription rates
row.translation=c() # Rows of output data that are translation rates
row.mRNA=c() # Rows of output data that are mRNA levels
row.protein=c() # Rows of output data that are protein levels
for(i in 1:ngene){
	row.transcription=c(row.transcription,2*i-1)
	row.translation=c(row.translation,2*i)
	row.mRNA=c(row.mRNA,2*i-1)
	row.protein=c(row.protein,2*i)
}

u=rep(1,2*ngene) # Mutation rates for each trait
sigma=rep(0.1,2*ngene) # SD of mutation effect size
width=rep(1,ngene) # Width of fitness function for each ngene
opt=rep(0,ngene)

T=1e4 # Duration of simulation for each lineage
Nrep=100 # Number of lineages

for(ncc in 1:length(coeff.all.m)){ # Go through the list of interaction matrices

coeff=coeff.all.m[[ncc]]
gt=list();pt=list()

for(n in 1:Nrep){
	gt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean genotypic values of the n-th lineage through time 
	pt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean phenotypes of the n-th lineage through time 
	pt[[n]][,1]=g2p(gt[[n]][,1],coeff) # Calculate the starting phenotype
	for(t in 2:(T+1)){
		gt[[n]][,t]=gt[[n]][,(t-1)]
		nm=rpois(1,lambda=sum(u)) # Total number of mutations that would occur in this time step
		if(nm>0){
			for(i in 1:nm){
				gt_ances=gt[[n]][,t]
				pt_ances=g2p(gt_ances,coeff) # Calcuate ancestral phenotype from the ancestral genotype
				gt_mutant=gt_ances
				type=sample(1:(2*ngene),1,prob=u/sum(u)) # Decide which trait the mutation affects
				effect=rnorm(1,mean=0,sd=sigma[type])
				gt_mutant[type]=gt_mutant[type]+effect # Add the mutation's effect to the affected trait
				pt_mutant=g2p(gt_mutant,coeff) # Calculate mutant phenotype from the mutant genotype
				pf=fix.prob(pt_ances[row.protein],pt_mutant[row.protein],width) # Fixation probability of the mutation
				if.fix=rbinom(n=1,size=1,prob=pf) # Decide if the mutation would fix
				if(if.fix==1){
					gt[[n]][,t]=gt_mutant # If the mutation is fixed, the mutant genotype would become the ancestral phenotype when the next mutation is examined
				}
			}
		}
		pt[[n]][,t]=g2p(gt[[n]][,t],coeff)
	}
}

gt_end=matrix(0,nrow=Nrep,ncol=6) # Genotypic values of all lineages at the end of simulation
pt_end=matrix(0,nrow=Nrep,ncol=6) # Phenotypes of all lineages at the end of simulation
# Create a new directory to store output files
dir_new=paste(ncc,"_",paste(width[1],width[2],width[3],sep=""),sep="")
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
	#setwd("..") # Back to parental directory before going to the next lineage
}
#setwd("..")
# Output end-point genotypic values and phenotypes
fn1=paste("gt_",paste(width[1],width[2],width[3],sep=""),".txt",sep="");write.table(gt_end,file=fn1,sep="\t")
fn2=paste("pt_",paste(width[1],width[2],width[3],sep=""),".txt",sep="");write.table(pt_end,file=fn2,sep="\t")
# Correlation matrix for all genotypic values and phenotypes
m.out=cov2cor(cov(data.frame(gt_end,pt_end)))
write.table(m.out,file="cor.txt",sep="\t")

# Data matrix that contains variances of all traits through time; the last column contains end-point variances that would be used in most analyses
var.all=matrix(0,nrow=ngene*4,ncol=(T+1))
for(t in 2:(T+1)){
	d=matrix(0,nrow=Nrep,ncol=ngene*4)
	for(n in 1:Nrep){
		d[n,1:6]=gt[[n]][,t]
		d[n,7:12]=pt[[n]][,t]
	}
	for(i in 1:(ngene*4)){
		var.all[i,t]=var(d[,i])
	}
}
write.table(var.all,file="var_all.txt",sep="\t")

}

