# Simulate coevolution of mRNA and protein levels
# Simplified model without considering protein degradation; degradation would be constant terms and are omitted
# Phenotypic values are all in log scale

Ne=1e3 # Effective population size
width=1 # SD of Gaussian fitness function (set width=0 for neutrality)

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

lambda.all=c(1,1) # Rates of mutations affecting transcription rate and (per-mRNA) translation rate
sig.all=(0.1,0.1) # Effect size SD of mutations affecting transcription rate and (per-mRNA) translation rate

T=1e4 # Duration of each simulation
Nrep=500 # Total number of simulations to run; can be interpreted as different genes or replicate lineages
v1=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time
v2=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time
v3=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time
out.end=matrix(0,nrow=Nrep,ncol=3) # All end-point values
for(n in 1:Nrep){
	v=matrix(0,nrow=3,ncol=(T+1)) # Phenotypes through time
	for(t in 2:(T+1)){
		v[,t]=v[,(t-1)]
		nm=rpois(1,lambda=sum(lambda.all))
		if(nm>0){
			for(i in 1:nm){
				phe_ances=v[3,t] # Ancestral protein level
				type=sample(1:2,1,prob=lambda.all/sum(lambda.all)) # Determine mutation type (which trait the mutation would affect)
				effect=rnorm(1,mean=0,sd=sig.all[type]) # Determine mutation effect size
				phe_mutant=phe_ances+effect # Mutant protein level
				pf=fix.prob(phe_ances,phe_mutant,width) # Fixation probability
				if.fix=rbinom(n=1,size=1,prob=pf) 
				v[type,t]=v[type,t]+effect*if.fix # Add the mutation's effect to the current transcription or translation rate if it is fixed
				v[3,t]=v[1,t]+v[2,t] # Update the protein level
			}
		}
		v1[n,t]=v[1,t];v2[n,t]=v[2,t];v1[n,t]=v[3,t] # Write the phenotype of lineage n at time t
	}
	out.end[n,]=v[,(T+1)] # Write end-point phenotype
}

# Variance through time (applicable if repeats are interpreted as lineages)
div=rep(0,nrow=3,ncol=(T+1))
for(t in 2:(T+1)){
	div[1,t]=var(v1[,t])
	div[2,t]=var(v2[,t])
	div[3,t]=var(v3[,t])
}

write.table(out.end,file="sim_basic_out_end.txt",sep="\t")
write.table(div,file="sim_basic_out_var.txt",sep="\t")

