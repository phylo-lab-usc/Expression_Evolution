#Modeling coevolution between mRNA and protein levels (simplified model without considering protein degradation)
#Two lineages simulated each time
#Protein level of the gene is under stabilizing selection in both lineages

Ne=1e4 #effective population size

# Fitness function
# Calculate fitness given distance to the optimal phenotype and shape of fitness function
fitness <- function(d,a){
	w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a) #Gaussian fitness function with SD equal to a
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

lambda1=1 # Rate of mutations affecting mRNA level
lambda2=1 # Rate of mutations affecting per-mRNA translation rate (translational efficiency)
sig1=0.1 # Variation of mutational effect on mRNA level (log scale)
sig2=0.1 # Variation of mutational effect on translational efficiency (log scale)
width=1 # Width of fitness fucntion (x-axis in log scale)

T=1e5 # Duration of the simulation

Ntest=500 # Number of independent simulations to run; each simulation can be interpreted as the representing a gene
# Matrices to store evolutionary divergence in mRNA level, translational efficiency and protein level, respectively
dist1=rep(0,Ntest);dist2=rep(0,Ntest);dist3=rep(0,Ntest)
for(test in 1:Ntest){
	v1=matrix(0,nrow=2,ncol=(T+1)) # Matrix to store mRNA levels of two lineages through time
	v2=matrix(0,nrow=2,ncol=(T+1)) # Matrix to store translational efficiency of two lineages through time
	v3=matrix(0,nrow=2,ncol=(T+1)) # Matrix to store protein levels of two lineages through time
	for(t in 2:(T+1)){
		for(n in 1:2){
			v1[n,t]=v1[n,(t-1)]
			v2[n,t]=v2[n,(t-1)]
			nm=rpois(1,lambda=(lambda1+lambda2)) # Total number of mutations that occurred in this time step
			if(nm>0){
				type=rbinom(n=nm,size=1,prob=lambda2/(lambda1+lambda2)) # Decide if the mutation affects mRNA level or translational efficiency
				for(i in 1:nm){
					phe_ances=v1[n,t]+v2[n,t] # Ancestral protein level
					if(type[i]==0){
						effect=rnorm(1,mean=0,sd=sig1) # Mutation's effect
						phe_mutant=phe_ances+effect # Mutant protein level
						pf=fix.prob(phe_ances,phe_mutant,width) # Calculate fixation probability
						if.fix=rbinom(n=1,size=1,prob=pf) # Decide whether the mutation would fix
						v1[n,t]=v1[n,t]+effect*if.fix # Add the mutation's effect to the population mean if it is fixed
					}else{
						effect=rnorm(1,mean=0,sd=sig2) # Mutation's effect
						phe_mutant=phe_ances+effect # Mutant protein level
						pf=fix.prob(phe_ances,phe_mutant,width) # Calculate fixation probability
						if.fix=rbinom(n=1,size=1,prob=pf) # Decide whether the mutation would fix
						v2[n,t]=v2[n,t]+effect*if.fix # Add the mutation's effect to the population mean if it is fixed
					}
				}
			}
			v3[n,t]=v1[n,t]+v2[n,t]
		}
	}
	#Calculate evolutionary divergence in each trait at the end of each simulation
	dist1[test]=v1[1,(T+1)]-v1[2,(T+1)]
	dist2[test]=v2[1,(T+1)]-v2[2,(T+1)]
	dist3[test]=v3[1,(T+1)]-v3[2,(T+1)]
}

#output
out=data.frame(dist1,dist2,dist3)
write.table(out,file="out.2sp.txt",sep="\t")

