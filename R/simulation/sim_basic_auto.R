# Modeling coevolution of mRNA and protein levels
# Simplified model without considering protein degradation; degradation would be constant terms and are omitted
# Phenotypic values are all in log scale
# Automatically simulate evolution under many combinations of effective population size and SD of fitness fucntion

# Parameter combinations to be considered
Ne.all=c(1e2,1e3,1e4,1e5) # Effective population size
width.all=c(1,5,10,50,100) # SD of Gaussian fitness function
par.all=matrix(0,nrow=length(Ne.all)*length(width.all),ncol=2)
row=0
for(i in 1:length(Ne.all)){
	for(j in 1:length(width.all)){
		row=row+1
		par.all[row,1]=Ne.all[i]
		par.all[row,2]=width.all[j]
	}
}

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

u.all=c(5e-4,5e-4) # Per-genome rates of mutations affecting transcription rate and (per-mRNA) translation rate
sig.all=c(0.1,0.1) # SD of effect size of mutations affecting transcription rate and (per-mRNA) translation rate

T=1e5 # Duration of each simulation; set to be equal to the biggest Ne
T.rec=(1:(T/100))*100
Nrep=500 # Total number of simulations to run; can be interpreted as different genes or replicate lineages

# Data matrices to store the output
out.end.all=matrix(0,nrow=nrow(par.all)*Nrep,ncol=5)
out.var.all=matrix(0,nrow=nrow(par.all)*T/100,ncol=6)

for(npar in 1:nrow(par.all)){
	Ne=par.all[npar,1]
	width=par.all[npar,2]
	lambda.all=u.all*2*Ne # Total rate of mutations (2*Ne*u)

	# Determine which rows in the output data matrix is for this parameter combination
	row.start=(npar-1)*Nrep+1
	row.end=npar*Nrep
	out.end.all[row.start:row.end,1]=Ne
	out.end.all[row.start:row.end,2]=width
	
	v1=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time
	v2=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time
	v3=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time
	for(n in 1:Nrep){
		vr=matrix(0,nrow=3,ncol=(T+1)) # Data matrix for this lineage's phenotype through time
		for(t in 2:(T+1)){
			vr[,t]=vr[,(t-1)]
			nm=rpois(1,lambda=sum(lambda.all))
			if(nm>0){
				for(i in 1:nm){
					phe_ances=vr[1,t]+vr[2,t] # Ancestral protein level
					type=sample(1:2,1,prob=lambda.all/sum(lambda.all)) # Determine whether the mutation affects transcription or translation
					effect=rnorm(1,mean=0,sd=sig.all[type]) # Determine mutation effect size
					phe_mutant=phe_ances+effect # Mutant protein level
					pf=fix.prob(phe_ances,phe_mutant,width,Ne) # Calculate fixation probability
					if.fix=rbinom(n=1,size=1,prob=pf) # Determine if the mutation gets fixed
					vr[type,t]=vr[type,t]+effect*if.fix # If the mutation if fixed, add its effect to the population mean phenotype
					vr[3,t]=vr[1,t]+vr[2,t] # Protein level after the mutation is considered
				}
			}
			v1[n,t]=vr[1,t]
			v2[n,t]=vr[2,t]
			v3[n,t]=vr[3,t]
		}
		out.end.all[(row.start+n-1),3:5]=vr[,(T+1)] # End-point phenotype for the replicate lineage
	}

	# Variance through time for the parameter value combination (write every 100 time steps)
	for(i in 1:(T/100)){
		t=i*100
		rowt=(npar-1)*(T/100)+i
		out.var.all[rowt,1]=Ne
		out.var.all[rowt,2]=width
		out.var.all[rowt,3]=t
		out.var.all[rowt,4]=var(v1[,(t+1)])
		out.var.all[rowt,5]=var(v2[,(t+1)])
		out.var.all[rowt,6]=var(v3[,(t+1)])
	}
}

# Output
write.table(out.end.all,file="out_basic_end_all.txt",sep="\t") # End-point phenotypes
write.table(out.var.all,file="out_basic_var_all.txt",sep="\t") # Variance through time
