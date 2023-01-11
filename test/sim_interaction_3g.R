# Modeling coevolution between mRNA and protein levels for three interacting genes
# Automatically run many manually chosen combinations of interaction parameters

setwd("/Users/daohanji/Desktop/Expression_Evolution")

dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
ngene=3 # Number of genes under consideration

# Interaction parameters
# Numbers represent elements [1,2], [1,3], [2,1], [2,3], [3,1], [3,2], respectively
par.all=rbind(c(0.5,0,0,0.5,0,0), # G1 activates G2, G2 activates G3
	c(0.5,0,0,0.5,-0.5,0), # G1 activates G2, G2 activates G3, G3 represses G1 (row 1 + negative feedback)
	c(0,0.5,0,0.5,0,0), # G1 and G2 independently activate G3
	c(0,0.5,0,0.5,-0.5,-0.5), # G1 and G2 independently activate G3 and both are repressed by G3 (row 3 + negative feedback)
	c(-0.5,0.5,-0.5,0.5,0,0)) # G1 and G2 independently activate G3 and repress each other (row 3 + "competition")

# Covert vector of interaction paramters to matrix format
v2m <- function(v){
	m=matrix(0,nrow=3,ncol=3)
	m[1,2:3]=v[1:2]
	m[2,1]=v[3];m[2,3]=v[4]
	m[3,1:2]=v[5:6]
	return(m)
}

# Fitness (multivariate Gaussian fitness function)
# Calculate the overall fitness given a set of traits' values (d, a vector) and SDs of their respective fitness functions (a, also a vector) 
fitness <- function(d,a){
	w=rep(0,length(d))
	for(i in 1:length(d)){
		if(a[i]!=0){
			w[i]=dnorm(d[i],mean=0,sd=a[i])/dnorm(0,mean=0,sd=a[i]) # Calculate a fitness value for each trait
		}
	}
	w=w[which(a!=0)] # Remove traits that do not affect fitness
	w=prod(w)^(1/length(w)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
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

# Genotype-phenotype map (multiple genes)
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
	if(det(mat)!=0){ # Check if Ax=B is solvable; if not, return nothing (the output will have length=0)
		pt=solve(mat,right) # Solve Ax=B equation system
	}else{
		pt=rep(1e4,4) # When Ax=b is not solvable, assign a phenotypic value that leads to zero fitness (when SD of the fitness function takes the highest value considered in the study)
	}
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
sig.all=rep(0.1,2*ngene) # SD of mutation effect size

Ne=1e3 # Effective population size
width=c(0,0,1) # Width of fitness function for each ngene

T=1e4 # Duration of simulation for each lineage
T.rec=(1:(T/10))*10 # Time points at which results would be written into the output matrix
Nrep=100 # Number of lineages

out.all=matrix(0,nrow=nrow(par.all)*T/10,ncol=28)

row=1
for(npar in 1:nrow(par.all)){
	coeff=v2m(par.all[npar,])
	out.all[row:(row+T/10-1),1]=par.all[npar,1]
	out.all[row:(row+T/10-1),2]=par.all[npar,2]
	out.all[row:(row+T/10-1),3]=par.all[npar,3]
	out.all[row:(row+T/10-1),4]=par.all[npar,4]
	out.all[row:(row+T/10-1),5]=par.all[npar,5]
	out.all[row:(row+T/10-1),6]=par.all[npar,6]
	gt=list();pt=list()
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
					pf=fix.prob(pt_ances[row.protein],pt_mutant[row.protein],width,Ne) # Fixation probability of the mutation
					if.fix=rbinom(n=1,size=1,prob=pf) # Decide if the mutation would fix
					if(if.fix==1){
						gt[[n]][,t]=gt_mutant # If the mutation is fixed, the mutant genotype would become the ancestral phenotype when the next mutation is examined
					}
				}
			}
			pt[[n]][,t]=g2p(gt[[n]][,t],coeff)
		}
	}
	# Go through simulated data and extract information of interest
	for(i in 1:length(T.rec)){
		t=T.rec[i]
		out.all[row+i-1,7]=t
		d=matrix(0,nrow=Nrep,ncol=12)
		for(n in 1:Nrep){
			d[n,1:6]=gt[[n]][,(t+1)]
			d[n,7:12]=pt[[n]][,(t+1)]
		}
		# Write variances in genotypic values and phenotypes
		for(column in 8:19){
			out.all[row+i-1,column]=var(d[,column-7])
		}
		# Write correlations between genotypic values and between phenotypes
		# Write ratio of variances of mRNA and protein levels
		# Only end point values for each interaction parameter combination are calculated and written
		if(t==T){
			out.all[row+i-1,20]=cor(d[,1],d[,2])
			out.all[row+i-1,21]=cor(d[,3],d[,4])
			out.all[row+i-1,22]=cor(d[,5],d[,6])
			out.all[row+i-1,23]=cor(d[,7],d[,8])
			out.all[row+i-1,24]=cor(d[,9],d[,10])
			out.all[row+i-1,25]=cor(d[,11],d[,12])
			out.all[row+i-1,26]=out.all[row+i-1,15]/out.all[row+i-1,14]
			out.all[row+i-1,27]=out.all[row+i-1,17]/out.all[row+i-1,16]
			out.all[row+i-1,28]=out.all[row+i-1,19]/out.all[row+i-1,18]
		}
	}
	row=row+T/10
}

# Generate output file
write.table(data.frame(out.all),file="out_interaction_3g.txt",sep="\t")

# Generate a concise output file that contains end point values only
out.end=out.all[which(out.all[,7]==T),]
out.end=data.frame(out.end[,1:6],out.end[,8:28])
write.table(out.end,file="out_interaction_3g_end.txt",sep="\t")

