# Modeling coevolution between mRNA and protein levels for three interacting genes

library(gtools)

Ne=1e3 # Effective population size
dcr=1;dcp=1 # Decay rates of mRNA and protein (assumed to be the same for all genes)
ngene=3 # Number of genes under consideration

# Interaction parameters
# All permutations of 0, a positive number and a negative number
# Numbers represent elements [1,2], [1,3], [2,1], [2,3], [3,1], [3,2], respectively
permut.all=permutations(3,6,c(-0.5,0,0.5),repeats.allowed=TRUE)

# Covert vector of interaction paramters to matrix format
v2m <- function(v){
	m=matrix(0,nrow=3,ncol=3)
	m[1,2:3]=v[1:2]
	m[2,1]=v[3];m[2,3]=v[4]
	m[3,1:2]=v[5:6]
	return(m)
}

# Fitness function
# Calculate fitness given distance to the optimal phenotype and shape of fitness function
fitness <- function(d,a){
	w=rep(0,length(d))
	for(i in 1:length(d)){
		if(a[i]!=0){
			w[i]=dnorm(d[i],mean=0,sd=a[i])/dnorm(0,mean=0,sd=a[i]) # Calculate a fitness value for each trait
		}else{
			w[i]=1
		}
	}
	w=prod(w)^(1/length(which(a!=0))) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
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

u=rep(1,2*ngene) # Mutation rates for each trait
sigma=rep(0.1,2*ngene) # SD of mutation effect size
width=c(0,0,1) # Width of fitness function for each ngene
opt=rep(0,ngene)

T=1e4 # Duration of simulation for each lineage
T.rec=(1:(T/10))*10 # Time points at which results would be written into the output matrix
Nrep=100 # Number of lineages

out.all=matrix(0,nrow=nrow(permut.all)*T/10,ncol=28)

row=1
for(ncc in 1:nrow(permut.all)){
	coeff=v2m(permut.all[ncc,])
	out.all[row:(row+T/10-1),1]=permut.all[ncc,1]
	out.all[row:(row+T/10-1),2]=permut.all[ncc,2]
	out.all[row:(row+T/10-1),3]=permut.all[ncc,3]
	out.all[row:(row+T/10-1),4]=permut.all[ncc,4]
	out.all[row:(row+T/10-1),5]=permut.all[ncc,5]
	out.all[row:(row+T/10-1),6]=permut.all[ncc,6]
	gt=list();pt=list()
	for(n in 1:Nrep){
		gt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean genotypic values of the n-th lineage through time 
		pt[[n]]=matrix(0,nrow=2*ngene,ncol=(T+1)) # Mean phenotypes of the n-th lineage through time 
		pt[[n]][,1]=g2p(gt[[n]][,1],coeff) # Calculate the starting phenotype
		for(t in 2:(T+1)){
			gt[[n]][,t]=gt[[n]][,(t-1)]
			nm=rpois(1,lambda=sum(u)) # Total number of mutations that would occur in this time step
			if(nm>0){
				for(i in 1:nm){
					gt_ances=gt[[n]][,t]
					pt_ances=g2p(gt_ances,coeff) # Calcuate ancestral phenotype from the ancestral genotype
					gt_mutant=gt_ances
					type=sample(1:(2*ngene),1,prob=u/sum(u)) # Decide which trait the mutation affects
					effect=rnorm(1,mean=0,sd=sigma[type])
					gt_mutant[type]=gt_mutant[type]+effect # Add the mutation's effect to the affected trait
					pt_mutant=g2p(gt_mutant,coeff) # Calculate mutant phenotype from the mutant genotype
					pf=fix.prob(pt_ances[row.protein],pt_mutant[row.protein],width) # Fixation probability of the mutation
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
out.end=data.frame(out.end[,1:6],out.end[8:28])
write.table(out.end,file="out_interaction_3g_end.txt",sep="\t")
