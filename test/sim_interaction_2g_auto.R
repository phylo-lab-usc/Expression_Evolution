# Modeling coevolution between mRNA and protein levels for two interacting genes, one of which is subject to selection for an optimal protein level
# Only one gene's protein level is directly under selection
# Automatically run many combinations of interaction parameters

setwd("/Users/daohanji/Desktop/Expression_Evolution")

dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)

# A matrix that contains all combinations of interaction parameters to be examined
coeff.all=c(-9:9)*0.1
par.all=matrix(0,nrow=length(coeff.all)^2,ncol=2)
row=0
for(i in 1:length(coeff.all)){
	for(j in 1:length(coeff.all)){
		row=row+1
		par.all[row,1]=coeff.all[i]
		par.all[row,2]=coeff.all[j]
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
	if(a==0){ # Neutrality
		p=1/(2*Ne)
	}else{
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
	}
	return(p)
}

# Genotype-phenotype map (2 genes)
# Calculate the phenotype (equilibrium mRNA and protein abundances) given genotypic values and interaction parameters
# Input and output are both in log scale
# dcr and dcp are constants are thus not among the parameters
g2p <- function(gt,coeff){
	a1=gt[1];b1=gt[2];a2=gt[3];b2=gt[4]
	right=c(log(dcp)-b1,log(dcr)-a2,log(dcp)-b2,log(dcr)-a1)
	mat=rbind(c(1,-1,0,0),c(0,coeff[1],-1,0),c(0,0,1,-1),c(-1,0,0,coeff[2]))
	if(det(mat)!=0){ # Check if Ax=B is solvable; if not, return nothing (the output will have length=0)
		pt=solve(mat,right) # Solve Ax=b equation system
	}else{
		pt=rep(1e4,4) # When Ax=b is not solvable, assign a phenotypic value that leads to zero fitness (when SD of the fitness function takes the highest value considered in the study)
	}
	return(pt) # Return log scale phenotypes
}

# Mutation rate of each trait (equivalent to 2*Ne*u)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
lambda.all=c(1,1,1,1)

# Variance of mutation effect size (log scale)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
sig.all=c(0.1,0.1,0.1,0.1)

Ne=1e3 # Effective population size
width=1 # Width of fitness fucntion (x-axis in log scale)
T=1e4 # Duration of the simulation
T.rec=(1:(T/10))*10 # Time points at which results would be written (write every 10 time steps)
Nrep=100 # Number of replicate lineages

# Matrix to store all the output
# Columns: interaction parameters, time, variances of genotypic values, variances of phenotypes, transcriptipn-translation correlations, RNA-protein correlations, protein/RNA variance ratios
# Row number equals cum*T/10 because simulation results would be written every 10 time steps
out.all=matrix(0,nrow=nrow(par.all)*T/10,ncol=17)

row=1
for(c in 1:nrow(par.all)){
	coeff=par.all[c,]
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
					pf=fix.prob(pt_ances[4],pt_mutant[4],width,Ne)
					if.fix=rbinom(n=1,size=1,prob=pf)
					if(if.fix==1){
						gt[[n]][,t]=gt_mutant # The mutant genotypic values becomes the ancestral genotypic values before the next mutation would be considered
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
			# Because the first column of the data matrix represents time 0, phenotypes of time t are in column (t+1)
			d[n,1:4]=gt[[n]][,(t+1)] # Extract genotypic values of lineage n at time t
			d[n,5:8]=pt[[n]][,(t+1)] # Extract phenotype of lineage n at time t
		}
				
		# Write variances in genotypic values and phenotypes
		# Because the first columns are for interaction parameters and time, the 4th to 11th columns are used to store the varainces
		# 4:7: genotypic values; 8-11: phenotypes
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

# Generate output file
write.table(data.frame(out.all),file="out_interaction_2g_pt1.txt",sep="\t")

# Generate a concise output file that contains end point values only
out.end=out.all[which(out.all[,3]==T),]
out.end=data.frame(out.end[,1:2],out.end[,4:17])
write.table(out.end,file="out_interaction_2g_pt1_end.txt",sep="\t")
