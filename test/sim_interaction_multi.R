# Modeling coevolution between mRNA and protein levels for two or more interacting genes

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
ngene=4 # Number of genes under consideration

# Interaction parameters
# Entry [i,j] represents effect of gene i on gene j
coeff=matrix(0,nrow=ngene,ncol=ngene)
# Generate random numbers as interaction parameters
# Can be assigned in different ways depending on the question to be addressed
# Need to check solvability before calculating the phenotype
for(i in 1:ngene){
	for(j in 1:ngene){
		if(i!=j){
			coeff[i,j]=runif(n=1,min=-1,max=1)
		}
	}
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
	if(det(mat!=0)){
		pt=solve(mat,right) # Solve Ax=B equation system
		return(pt) # Return log scale phenotypes
	}
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

lambda.all=rep(1,2*ngene) # Mutation rates for each trait
sig.all=rep(0.1,2*ngene) # SD of mutation effect size for each trait

width=rep(1,ngene) # Width of fitness function for each gene
opt=matrix(0,nrow=n,ncol=ngene)

T=1e4 # Duration of simulation for each lineage
Nrep=100 # Number of lineages
gt=list() # Genotypic values of all lineages through time (each lineage would be a matrix in the list)
pt=list() # Phenotypes of all lineages through time (each lineage would be a matrix in the list)
gt_end=matrix(0,nrow=Nrep,ncol=2*ngene) # Genotypic values of all lineages at the end of simulation
pt_end=matrix(0,nrow=Nrep,ncol=2*ngene) # Phenotypes of all lineages at the end of simulation			
for(n in 1:Nrep){
	gt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean genotypic values of the n-th lineage through time 
	pt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean phenotypes of the n-th lineage through time 
	pt[[n]][,1]=g2p(gt[[n]][,1],coeff) # Calculate the starting phenotype
	for(t in 2:(T+1)){
		gt[[n]][,t]=gt[[n]][,(t-1)]
		nm=rpois(1,lambda=sum(lambda.all)) # Total number of mutations that would occur in this time step
		if(nm>0){
			for(i in 1:nm){
				gt_ances=gt[[n]][,t]
				pt_ances=g2p(gt_ances,coeff) # Calcuate ancestral phenotype from the ancestral genotype
				gt_mutant=gt_ances
				type=sample(1:(2*ngene),1,prob=lambda.all/sum(lambda.all)) # Decide which trait the mutation affects
				effect=rnorm(1,mean=0,sd=sig.all[type])
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
	gt_end[n,]=gt[[n]][,(T+1)]
	pt_end[n,]=pt[[n]][,(T+1)]
}

