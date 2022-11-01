# Modeling coevolution between mRNA and protein levels for two interacting genes, one of which is subject to selection for an optimal protein level
# Automatically run many combinations of interaction parameters

setwd("/Users/daohanji/Desktop/Expression_Evolution")

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)

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
T.rec=(1:(T/10))*10
Nrep=100 # Number of replicate lineages

coeff.all=c(-9:9)*0.1 # All interaction parameter values to be used; all combinations to be examined

# Number of interaction parameter combinations with redundant ones excluded
cnum=length(coeff.all)^2

# Matrix to store all the output
# Columns: interaction parameters, time, variances of genotypic values, variances of phenotypes, transcriptipn-translation correlations, RNA-protein correlations, protein/RNA variance ratios
# Row number equals cum*T/10 because simulation results would be written every 10 time steps
out.all=matrix(0,nrow=cnum*T/10,ncol=17)

row=1
for(c1 in 1:length(coeff.all)){ # effect of gene 1 on gene 2
	for(c2 in 1:length(coeff.all)){ # effect of gene 2 on gene 1
		coeff=c(coeff.all[c1],coeff.all[c2]) # Vector of interaction parameters
		out.all[row:(row+T/10-1),1]=coeff[1]
		out.all[row:(row+T/10-1),2]=coeff[2]
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
						effect=rnorm(1,mean=0,sd=sig1)
						gt_mutant[type]=gt_mutant[type]+effect
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
		for(i in 1:length(T.rec)){
			t=T.rec[i]
			out.all[row+i-1,3]=t
			d=matrix(0,nrow=Nrep,ncol=8)
			for(n in 1:Nrep){
				d[n,1:4]=gt[[n]][,(t+1)]
				d[n,5:8]=pt[[n]][,(t+1)]
			}
				
			# Write variances in genotypic values and phenotypes
			for(column in 4:11){
				out.all[row+i-1,column]=var(d[,column-3])
			}
				
			# Write correlations between genotypic values and between phenotypes
			# Write ratio of variances of mRNA and protein levels
			# Only end point values for each interaction parameter combination are calculated and written
			if(t==T){
				out.all[row+i-1,12]=cor(d[,1],d[,2])
				out.all[row+i-1,13]=cor(d[,3],d[,4])
				out.all[row+i-1,14]=cor(d[,5],d[,6])
				out.all[row+i-1,15]=cor(d[,7],d[,8])
				out.all[row+i-1,16]=out.all[row+i-1,9]/out.all[row+i-1,8]
				out.all[row+i-1,17]=out.all[row+i-1,11]/out.all[row+i-1,10]
			}
		}
		row=row+T/10
	}
}

# Generate output file
write.table(data.frame(out.all),file="out_interaction_2g_pt1.txt",sep="\t")

# Generate a concise output file that contains end point values only
out.end=out.all[which(out.all[,3]==T),]
out.end=data.frame(out.end[,1:2],out.end[,4:17])
write.table(out.end,file="out_interaction_2g_pt1_end.txt",sep="\t")
