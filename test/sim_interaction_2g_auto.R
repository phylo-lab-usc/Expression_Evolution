# Modeling coevolution between mRNA and protein levels for two interacting genes, one of which is subject to selection for an optimal protein level
# Only one gene's protein level is directly under selection
# Automatically run many combinations of interaction parameters

dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
coeff.all=c(-9:9)*0.1 # All interaction parameter values to be used; all combinations to be examined
cnum=length(coeff.all)^2 # Number of interaction parameter combinations (no redundant set as only one gene is under direct selection)

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
fix.prob <- function(x1,x2,a,Ne){
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
# Calculate the phenotype (equilibrium mRNA and protein abundances) given genotypic values and interaction parameters
# Input and output are both in log scale
g2p <- function(gt,coeff){
	a1=gt[1];b1=gt[2];a2=gt[3];b2=gt[4]
	right=c(log(dcp)-b1,log(dcr)-a2,log(dcp)-b2,log(dcr)-a1)
	mat=rbind(c(1,-1,0,0),c(0,coeff[1],-1,0),c(0,0,1,-1),c(-1,0,0,coeff[2]))
	pt=solve(mat,right)
	return(pt) # Return log scale phenotypes
}

# Mutation rate of each trait (Ne*u)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
lambda.all=c(1,1,1,1)

# Variance of mutation effect size (log scale)
# For transcription rate of gene 1, translation rate of gene 1, transcription rate of gene 2, translation rate of gene 2, respectively 
sig.all=c(0.1,0.1,0.1,0.1)

Ne=1e3 # Effective population size
width=1 # Width of fitness fucntion (x-axis in log scale)
T=1e4 # Duration of the simulation
T.rec=(1:(T/10))*10 # Time points at which results would be written
Nrep=100 # Number of replicate lineages

# Matrix to store all the output
# Columns: interaction parameters, time, variances of genotypic values, variances of phenotypes, transcriptipn-translation correlations, RNA-protein correlations, protein/RNA variance ratios
# Row number equals cum*T/10 because simulation results would be written every 10 time steps
out.all=matrix(0,nrow=cnum*T/10,ncol=17)

row=1
for(c1 in 1:length(coeff.all)){ # Effect of gene 1 on gene 2
	for(c2 in 1:length(coeff.all)){ # Effect of gene 2 on gene 1
		coeff=c(coeff.all[c1],coeff.all[c2]) # Vector of interaction parameters
		# Write interaction parameters
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
						type=sample(1:4,1,prob=lambda.all/sum(lambda.all)))
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
