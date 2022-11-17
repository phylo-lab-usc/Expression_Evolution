# Modeling coevolution between mRNA and protein levels for two interacting genes that are both under selection
# Both genes are under direct selection
# Automatically run many combinations of interaction parameters

setwd("/Users/daohanji/Desktop/Expression_Evolution")

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
ngene=2 # Number of genes under consideration
# All values of the interaction parameter
coeff.all=c(-9:9)*0.1
# Number of interaction parameter combinations with redundant ones excluded
cnum=length(coeff.all)^2/2+length(coeff.all)/2

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
	w=prod(w)^(1/length(w)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits (neutral traits excluded)
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
	pt=solve(mat,right) # Solve Ax=B equation system
	return(pt) # Return log scale phenotypes
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

T=1e4 # Duration of simulation for each lineage
T.rec=(1:(T/10))*10 # Time points at which results would be written into the output matrix
Nrep=100 # Number of lineages

# Matrix to store all the output
# Columns: interaction parameters, time, variances of genotypic values, variances of phenotypes, transcriptipn-translation correlations, RNA-protein correlations, protein/RNA variance ratios
# Row number equals cum*T/10 because simulation results would be written every 10 time steps
out.all=matrix(0,nrow=cnum*T/10,ncol=17)

row=1
for(c1 in 1:length(coeff.all)){
	for(c2 in 1:length(coeff.all)){
		if(c1<=c2){ # Remove redundant combinations
			# Generate the interaction matrix
			coeff=matrix(0,nrow=2,ncol=2)
			coeff[1,2]=coeff.all[c1]
			coeff[2,1]=coeff.all[c2]
			# Write interaction parameters
			out.all[row:(row+T/10-1),1]=coeff[1,2]
			out.all[row:(row+T/10-1),2]=coeff[2,1]
			gt=list() # Genotypic values of all lineages through time (each lineage would be a matrix in the list)
			pt=list() # Phenotypes of all lineages through time (each lineage would be a matrix in the list)
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
					pt[[n]][,t]=g2p(gt[[n]][,t],coeff) # Calculate population mean phenotypes from genotypic values
				}
			}

			# Go through simulated data and extract information of interest
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
}
# Generate output file
write.table(data.frame(out.all),file="out_interaction_2g_pt2.txt",sep="\t")

# Generate a concise output file that contains end point values only
out.end=out.all[which(out.all[,3]==T),]
out.end=data.frame(out.end[,1:2],out.end[4:17])
write.table(out.end,file="out_interaction_2g_pt2_end.txt",sep="\t")



