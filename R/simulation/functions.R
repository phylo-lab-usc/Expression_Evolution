# Functions that would be used in most simulations

# Fitness (univariate Guassion fitness function)
# Calculate fitness given distance to optimum (d) and SD of the Gaussian fitness function (a)
# The result remains the same if d and a are multiplied by the same number
fitness <- function(d,a){
	if(a==0){ # Neutrality
		w=1
	}else{
		w=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a)
	}
	return(w)
}

# An alternative way to write the above function
fitness <- function (d,a){
	if(a==0){ # Neutrality
		w=1
	}else{
		w=exp(-0.5*(d^2/a^2))
	}
	return(w)
}

# Fitness (multivariate Gaussian fitness function)
# Calculate the overall fitness given a set of traits' values (d, a vector) and SDs of their respective fitness functions (a, also a vector) 
# Overall fitness calculated as geometric mean of fitness values calculated from individual traits
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

# Fitness (multivariate Gaussian fitness function)
# Calculate the overall fitness given a set of traits' values (d, a vector) and SDs of their respective fitness functions (a, also a vector) 
# Overall fitness calculated from Euclidean distance from the global optimum; fitness function's SD becomes a scaling coefficient (remains mathematically equivalent to original definition)
fitness <- function(d,a){
	# Remove traits that do not affect fitness
	ds=d[which(a>0)]
	as=a[which(a>0)]
	ds=ds/as # Rescale the phenotypes; when the fitness function is narrow, the distance from optimum is rescaled to be larger
	D=sqrt(sum(ds^2)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
	w=dnorm(D,mean=0,sd=1)/dnorm(0,mean=0,sd=1)
	return(w)
}

# Fitness with cost of expression taken into account (univariate)
# Parameters: the phenotype (mRNA level and the protein level, log scale), the optimum, per-molecule cost of the mRNA, per-molecule cost of the protein
fitness <- function(v,opt,a,cr,cp){
	# Fitness with respect to the protein level (cost not considered)
	if(a==0){ # Neutrality
		wp=1
	}else{
		d=v[2]-opt
		wp=dnorm(d,mean=0,sd=a)/dnorm(0,mean=0,sd=a)
	}
	cost=cr*exp(v[1])+cp*exp(v[2]) # Fitness reduction due to cost of expression
	w=wp-cost # Realized fitness 
	if(w<0){
		w=0 # Set fitness to 0 if production cost exceeds benefit of expression
	}
	return(w)
}

# Fitness with cost of expression taken into account (multivariate)
# Parameters: the phenotype (mRNA level and the protein level, log scale), the optimum, per-molecule cost of the mRNA, per-molecule cost of the protein
# vr, vp, opt and a are all vectors
# Cost is assumed to be equal for all genes, so cr and cp are constant
fitness <- function(vr,vp,opt,a,cr,cp){
	# Fitness with respect to the protein levels (cost not considered)
	d=vp-opt
	ds=d[which(a>0)]
	as=a[which(a>0)]
	ds=ds/as # Rescale the phenotypes; when the fitness function is narrow, the distance from optimum is rescaled to be larger
	D=sqrt(sum(ds^2)) # Overall fitness calculated as geometric mean of fitness values calculated from individual traits
	wp=dnorm(D,mean=0,sd=1)/dnorm(0,mean=0,sd=1)

	cost=sum(cr*exp(vr))+sum(cp*exp(vp)) # Fitness reduction due to cost of expression
	w=wp-cost # Realized fitness 
	if(w<0){
		w=0 # Set fitness to 0 if production cost exceeds benefit of expression
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

