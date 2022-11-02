# Modeling coevolution of mRNA and protein levels
# Simplified model without considering protein degradation; degradation would be constant terms and are omitted
# Phenotypic values are all in log scale

Ne=1e3 # Effective population size
width=1 # SD of Gaussian fitness function

#fitness function
fitness <- function(d,a){ #parameters: distance to optimum and width of fitness function
	w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a)
	return(w)
}

#fixation probability
fix.prob <- function(x1,x2,a){
	wa=fitness(x1,a)
	wm=fitness(x2,a)
	if(wa>0){
		s=(wm/wa)-1
		if(s==0){
			p=1/(2*Ne)
		}else{
			p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
		}
	}else{
		if(wm==0){ #both ancestral and mutant fitness are close to 0
			p=1/(2*Ne)
		}else{ #ancestral fitness is close to 0 while mutant fitness isn't
			p=1
		}
	}
	return(p)
}

lambda1=1 # Rate of mutations affecting transcription rate
lambda2=1 # Rate of mutations affecting (per-mRNA) translation rate
sig1=0.1 # SD of effect size of mutations affecting transcription rate
sig2=0.1 # SD of effect size of mutations affecting (per-mRNA) translation rate

T=1e5 # Duration of each simulation
Nrep=500 # Total number of simulations to run; can be interpreted as different genes or replicate lineages
v1=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time
v2=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time
v3=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time
out.end=matrix(0,nrow=Nrep,ncol=3) # All end-point values
for(n in 1:Nrep){
	for(t in 2:(T+1)){
		v1[n,t]=v1[n,(t-1)]
		v2[n,t]=v2[n,(t-1)]
		nm=rpois(1,lambda=(lambda1+lambda2))
		if(nm>0){
			type=rbinom(n=nm,size=1,prob=lambda2/(lambda1+lambda2))
			for(i in 1:nm){
				phe_ances=v1[n,t]+v2[n,t]
				if(type[i]==0){
					effect=rnorm(1,mean=0,sd=sig1)
					phe_mutant=phe_ances+effect
					pf=fix.prob(phe_ances,phe_mutant,width)
					if.fix=rbinom(n=1,size=1,prob=pf)
					v1[n,t]=v1[n,t]+effect*if.fix
				}else{
					effect=rnorm(1,mean=0,sd=sig2)
					phe_mutant=phe_ances+effect
					pf=fix.prob(phe_ances,phe_mutant,width)
					if.fix=rbinom(n=1,size=1,prob=pf)
					v2[n,t]=v2[n,t]+effect*if.fix
				}
			}
		}
		v3[n,t]=v1[n,t]+v2[n,t]
	}
	out.end[n,1]=v1[n,(T+1)]
	out.end[n,2]=v2[n,(T+1)]
	out.end[n,3]=v3[n,(T+1)]
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



