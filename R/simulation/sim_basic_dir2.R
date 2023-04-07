# Modeling coevolution of mRNA and protein levels
# Simplified model without considering protein degradation; degradation would be constant terms and are omitted
# Phenotypic values are all in log scale
# Simulate multiple genes with different optimal protein levels

# Simulate evolution of multiple genes with different optimal protein levels

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

ngene=20 # Number of genes to simulate
Ne=1e3 # Effective population size
width=1 # SD of Gaussian fitness function (assumed to be constant for different genes)
lambda.allgene=matrix(1,nrow=ngene,ncol=2) # Rates of mutations affecting transcription rate and (per-mRNA) translation rate (equivalent to Ne*u)
sig.allgene=matrix(0.1,nrow=ngene,ncol=2) # SD of effect size of mutations affecting transcription rate and (per-mRNA) translation rate
opt.all=rnorm(ngene,mean=0,sd=2) # Optimal protein levels of the genes (can be specified in different ways)

T=1e5 # Duration of each simulation
T.rec=(1:(T/100))*100
Nrep=500 # Total number of simulations to run; can be interpreted as different genes or replicate lineages
out.end.all=matrix(0,nrow=Nrep*ngene,ncol=4) # All end-point values
var.all=matrix(0,nrow=length(T.rec)*Nrep,ncol=5)
for(g in 1:ngene){
	lambda.all=lambda.allgene[g,]
	sig.all=sig.allgene[g,]
	opt=opt.all[g]
	row.start=(g-1)*Nrep+1
	row.end=g*Nrep
	out.end.all[row.start:row.end,1]=opt
	
	v1=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for transcription rates through time
	v2=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for translation rates through time
	v3=matrix(0,nrow=Nrep,ncol=(T+1)) # Data matrix for protein levels through time
	out.end=matrix(0,nrow=Nrep,ncol=3) # All end-point values
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
					pf=fix.prob(phe_ances-opt,phe_mutant-opt,width,Ne) # Calculate fixation probability
					if.fix=rbinom(n=1,size=1,prob=pf) # Determine if the mutation gets fixed
					vr[type,t]=vr[type,t]+effect*if.fix # If the mutation if fixed, add its effect to the population mean phenotype
					vr[3,t]=vr[1,t]+vr[2,t] # Protein level after the mutation is considered
				}
			}
			v1[n,t]=vr[1,t]
			v2[n,t]=vr[2,t]
			v3[n,t]=vr[3,t]
		}
		out.end[n,]=vr[,(T+1)] # End-point phenotype for the replicate lineage
	}
	out.end.all[row.start:row.end,2:4]=out.end

	# Variance through time (applicable if repeats are interpreted as lineages)
	var.row.start=(g-1)*length(T.rec)+1
	var.row.end=g*length(T.rec)
	var.all[var.row.start:var.row.end,1]=opt
	div=matrix(0,nrow=length(T.rec),ncol=4)
	for(i in 1:length(T.rec)){
		t=T.rec[i]
		div[i,1]=t
		div[i,2]=var(v1[,t])
		div[i,3]=var(v2[,t])
		div[i,4]=var(v3[,t])
	}
	var.all[var.row.start:var.row.end,2:5]=div
}

# Output
# Fitness function SD is included in the file name
fn1=paste("multi_",width,"_",T,"_out_basic_end.txt",sep="") # End-point phenotype
write.table(out.end.all,file=fn1,sep="\t")
fn2=paste("multi_",width,"_",T,"_out_basic_var.txt",sep="") # Variance through time
write.table(var.all,file=fn2,sep="\t")



















