setwd("/Users/daohanji/Desktop/Expression_Evolution/out_basic")

d<-read.table("0_out_basic_end.txt",sep="\t")
#de=d[which(d[,1]!=0&d[,2]!=0),] # Remove lineages that with zero evolutionary change in either trait

plot(d[,1],d[,2],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
abline(lm(d[,2]~d[,1]),col="blue")

plot(d[,1],d[,3],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
abline(lm(d[,3]~d[,1]),col="blue")

