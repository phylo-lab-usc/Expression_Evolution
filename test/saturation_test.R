# Check whether variance of a trait under selection is saturated at the end of the simulation
# For single case examination, not for mass test

width=1
coeff=c(0,0)

dir=paste(width,"_",coeff[1],"_",coeff[2])
setwd(dir)

fn=paste("var_all_",width,"_",coeff[1],"_",coeff[2],".txt",sep="")
v <- read.table(fn,sep="\t")

# Decide which row to examine
# Should be 8 for the case of two genes with one affecting fitness
row.test=8
# Calculate rate of change
dv=c()
for(t in 2:ncol(v)){
	dv=c(dv,v[row.test,t]-v[row.test,(t-1)])
}

# Examine distribution of per-generation change in variance for different time intervals (custom defined)
wilcox.test(dv[1:5000],mu=0) # Test whether change of variance per-generation is significantly different from 0 for the first half of time
wilcox.test(dv[5001:10000],mu=0) # Test whether change of variance per-generation is significantly different from 0 for the later half of time

