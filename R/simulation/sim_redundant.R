# Modeling coevolution of mRNA and protein levels
# Simplified model without considering protein degradation; degradation would be constant terms and are omitted
# Phenotypic values are all in log scale

# Simulating two functionally equivalent genes (fitness determined by sum of their protein levels)

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

Ne=1e3 # Effective population size
width=1 # SD of Gaussian fitness function
lambda.all=c(1,1,1,1) # Rates of mutations affecting transcription rate and (per-mRNA) translation rate
sig.all=c(0.1,0.1,0.1,0.1) # SD of effect size of mutations affecting transcription rate and (per-mRNA) translation rate

T=1e5 # Duration of each simulation
Nrep=500 # Total number of simulations to run; can be interpreted as different genes or replicate lineages
v1=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time (gene 1)
v2=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time (gene 1)
v3=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time (gene 2)
v4=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time (gene 2)
v5=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time (gene 1)
v6=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time (gene 2)
out.end=matrix(0,nrow=Nrep,ncol=6) # All end-point values
for(n in 1:Nrep){
	vr=matrix(0,nrow=4,ncol=(T+1)) # Data matrix for this lineage's transcription and translation rates through time
	for(t in 2:(T+1)){
		vr[,t]=vr[,(t-1)]
		nm=rpois(1,lambda=sum(lambda.all))
		if(nm>0){
			for(i in 1:nm){
				phe_ances=vr[,t]
				total_ances=log(exp(phe_ances[1]+phe_ances[2])+exp(phe_ances[3]+phe_ances[4])) # Total protein level of the ancestor
				type=sample(1:4,1,prob=lambda.all/sum(lambda.all)) # Determine which trait is affected by the mutation
				effect=rnorm(1,mean=0,sd=sig.all[type]) # Determine mutation effect size
				phe_mutant=phe_ances;phe_mutant[type]=phe_mutant[type]+effect
				total_mutant=log(exp(phe_mutant[1]+phe_mutant[2])+exp(phe_mutant[3]+phe_mutant[4])) # Total protein level of the mutant
				pf=fix.prob(total_ances,total_mutant,width,Ne) # Determine if the mutation gets fixed
				if.fix=rbinom(n=1,size=1,prob=pf) # If the mutation if fixed, add its effect to the population mean phenotype
				vr[type,t]=vr[type,t]+effect*if.fix
			}
		}
		v1[n,t]=vr[1,t]
		v2[n,t]=vr[2,t]
		v3[n,t]=vr[3,t]
		v4[n,t]=vr[4,t]
		v5[n,t]=vr[1,t]+vr[2,t]
		v6[n,t]=vr[3,t]+vr[4,t]
	}
	out.end[n,1:4]=vr[,(T+1)] # End-point phenotype for the replicate lineage
}
out.end[,5]=out.end[,1]+out.end[,2]
out.end[,6]=out.end[,3]+out.end[,4]

# Variance through time (applicable if repeats are interpreted as lineages)
div=matrix(0,nrow=6,ncol=(T+1))
for(t in 2:(T+1)){
	div[1,t]=var(v1[,t])
	div[2,t]=var(v2[,t])
	div[3,t]=var(v3[,t])
	div[4,t]=var(v4[,t])
	div[5,t]=var(v5[,t])
	div[6,t]=var(v6[,t])
}

# Output
# Fitness function SD is included in the file name (0 corresponds to neutrality)
fn1=paste(width,"_out_redundant_end.txt",sep="") # End-point phenotype of each replicate lineage
write.table(out.end,file=fn1,sep="\t")
fn2=paste(width,"_out_redundant_var.txt",sep="") # Variance through time
write.table(div,file=fn2,sep="\t")
