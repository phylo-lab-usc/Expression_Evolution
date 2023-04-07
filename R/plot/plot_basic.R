# Plotting results for evolution on replicate lineages with stabilizing selection on the protein level or under neutrality

library(ggplot2)

# Time-variance plot (stabilizing selection)
d<-read.table("out_basic_var_all.txt",sep="\t") # Data file: variance through time under stabilizing selection
dm=d[which(d[,1]==1e3&d[,2]==1),] # Keep rows corresponding to "default" values of Ne and fitness function SD

# Reorganize the data to prepare them for ggplot2
col1=rep(dm[,3],3) # Column of time
col2=c(rep("mRNA",nrow(dm)),rep("Translation",nrow(dm)),rep("Protein",nrow(dm))) # Column of trait
col3=c(dm[,4],dm[,5],dm[,6]) # Column of variance
dnew=data.frame(col1,col2,col3)
colnames(dnew)=c("t","Trait","div")
dnew$Trait=factor(dnew$Trait,levels=c("mRNA","Translation","Protein"),ordered=TRUE) # Order trait names to be shown in the color legend
g<-ggplot(dnew,aes(x=t,y=div,colour=Trait),pallete=c())+geom_point() # Make scatterplot
g=g+theme_classic() # Clear background
g=g+xlab("Divergence time")+ylab("Variance") # Axis labels
g=g+scale_color_manual(values=c("orange","grey","purple")) # Color data points for different traits differently
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font size of axis labels, axis text, legend title, and legend text
ggsave("t-var.pdf",plot=g,width=7.5,height=5) # Save the plot as a file

# Use basic plot function
#plot(dm[,3],dm[,4],xlab="T",ylab="Transcription rate (mRNA level) variance")
#plot(dm[,3],dm[,5],xlab="T",ylab="Translation rate variance")
#plot(dm[,3],dm[,6],xlab="T",ylab="Protein level variance")
#plot(dm[1:200,3],dm[1:200,6],xlab="T",ylab="Protein level variance")
#plot(dm[,3],dm[,6]/dm[,4],xlab="T",ylab="Variance ratio")

# Time-variance plot (neutral)
dm<-read.table("0_out_basic_var.txt",sep="\t") # Data file: variance through time under neutrality

# Reorganize the data to prepare them for ggplot2
col1=rep(dm[,1],1) # Column of time
col2=c(rep("mRNA",nrow(dm)),rep("Translation",nrow(dm)),rep("Protein",nrow(dm))) # Column of trait
col3=c(dm[,2],dm[,3],dm[,4]) # Column of variance
dnew=data.frame(col1,col2,col3)
colnames(dnew)=c("t","Trait","div")
dnew$Trait=factor(dnew$Trait,levels=c("mRNA","Translation","Protein"),ordered=TRUE) # Order trait names to be shown in the color legend
g<-ggplot(dnew,aes(x=t,y=div,colour=Trait),pallete=c())+geom_point() # Make scatterplot
g=g+theme_classic() # Clear background
g=g+xlab("Divergence time")+ylab("Variance") # Axis labels
g=g+scale_color_manual(values=c("orange","grey","purple")) # Color data points for different traits differently
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font size of axis labels, axis text, legend title, and legend text
ggsave("t-var_neutral.pdf",plot=g,width=7.5,height=5)  # Save the plot as a file

# Use basic plot function
#plot(dm[,1],dm[,2],xlab="T",ylab="Transcription rate (mRNA level) variance")
#plot(dm[,1],dm[,3],xlab="T",ylab="Translation rate variance")
#plot(dm[,1],dm[,4],xlab="T",ylab="Protein level variance")
#plot(dm[,1],dm[,4]/dm[,2],xlab="T",ylab="Variance ratio")

# End-point result
d<-read.table("out_basic_end_all.txt",sep="\t") # Data file: end-point phenotypes of replicate lineages that evolved under stabilizing selection
dm=d[which(d[,1]==1e3&d[,2]==1),] # Keep rows corresponding to "default" values of Ne and fitness function SD

dnew=data.frame(dm[,3],dm[,4]) # Extract traits of interest
colnames(dnew)=c("mRNA","translation")
g<-ggplot(dnew,aes(x=mRNA,y=translation))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Translation rate") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor1-basic.pdf",plot=g,width=5.5,height=5)  # Save the plot as a file

dnew=data.frame(dm[,3],dm[,5]) # Extract traits of interest
colnames(dnew)=c("mRNA","protein")
g<-ggplot(dnew,aes(x=mRNA,y=protein))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Protein level") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor2-basic.pdf",plot=g,width=5.5,height=5)  # Save the plot as a file

# Use basic plot function
#plot(dm[,3],dm[,4],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
#abline(lm(dm[,4]~dm[,3]),col="blue",lwd=2)

#plot(dm[,3],dm[,5],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
#abline(lm(dm[,5]~dm[,3]),col="blue")

# End-point result (neutral)
dm<-read.table("0_out_basic_end.txt",sep="\t") # Data file: end-point phenotypes of replicate lineages that evolved under neutrality

dnew=data.frame(dm[,1],dm[,2]) # Extract traits of interest
colnames(dnew)=c("mRNA","translation")
g<-ggplot(dnew,aes(x=mRNA,y=translation))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Translation rate") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor1-basic-neutral.pdf",plot=g,width=5.5,height=5)  # Save the plot as a file

dnew=data.frame(dm[,1],dm[,3]) # Extract traits of interest
colnames(dnew)=c("mRNA","protein")
g<-ggplot(dnew,aes(x=mRNA,y=protein))
g=g+geom_point()+geom_smooth(method="lm") # Make scatterplot and add least-squares regression line
g=g+theme_classic() # Clear background
g=g+xlab("Transcription rate (mRNA level)")+ylab("Protein level") # Axis labels
g=g+theme(axis.text=element_text(size=12),axis.title=element_text(size=15)) # Adjust font size of axis labels and text
ggsave("cor2-basic-neutral.pdf",plot=g,width=5.5,height=5)  # Save the plot as a file

# Use basic plot function
#plot(dm[,1],dm[,2],xlab="Transcription rate (mRNA level)",ylab="Translation rate",cex.lab=1.5)
#abline(lm(dm[,2]~dm[,1]),col="blue",lwd=2)

#plot(dm[,1],dm[,3],xlab="Transcription rate (mRNA level)",ylab="Protein level",cex.lab=1.5)
#abline(lm(dm[,3]~dm[,1]),col="blue")

