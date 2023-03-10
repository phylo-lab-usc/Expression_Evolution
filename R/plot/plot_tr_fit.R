#setwd("/Users/daohanji/Desktop/Expression_Evolution/rerun/out_tree/error")
setwd("/Users/rexjiang/Desktop/Lab/multivariate_trait/expression/rerun/out_tree")

library(ggplot2)

# Stabilizing selection
d<-read.table("1_out_sum_error_tr_new.txt",sep="\t")

# Reorganize the data to prepare them for ggplot2
col1=rep(d[,1],2)
col2=c(rep("BM",6),rep("OU",6))
dm=data.frame(col1,col2,c(d[,2],d[,3]),c(d[,4],d[,5]),c(d[,6],d[,7]))
colnames(dm)=c("esd","model","w1","w2","w3")

g1<-ggplot(dm,aes(x=esd,y=w1,fill=model))+geom_bar(stat="identity")
g1=g1+theme_classic() # Clear background
g1=g1+scale_fill_manual(values=c("orange","purple")) # Set the colors
g1=g1+xlab("Error SD")+ylab("AIC weight") # Axis labels
g1=g1+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_mRNA.pdf",plot=g1,width=8,height=5) # Save the plot

g2<-ggplot(dm,aes(x=esd,y=w2,fill=model))+geom_bar(stat="identity")
g2=g2+theme_classic() # Clear background
g2=g2+scale_fill_manual(values=c("orange","purple")) # Set the colors
g2=g2+xlab("Error SD")+ylab("AIC weight") # Axis labels
g2=g2+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_translation.pdf",plot=g2,width=8,height=5) # Save the plot

g3<-ggplot(dm,aes(x=esd,y=w3,fill=model))+geom_bar(stat="identity")
g3=g3+theme_classic() # Clear background
g3=g3+scale_fill_manual(values=c("orange","purple")) # Set the colors
g3=g3+xlab("Error SD")+ylab("AIC weight") # Axis labels
g3=g3+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_protein.pdf",plot=g3,width=8,height=5) # Save the plot

# Results under neutrality
# Reorganize the data to prepare them for ggplot2
d<-read.table("0_out_sum_error_tr_new.txt",sep="\t")
col1=rep(d[,1],2)
col2=c(rep("BM",6),rep("OU",6))
dm=data.frame(col1,col2,c(d[,2],d[,3]),c(d[,4],d[,5]),c(d[,6],d[,7]))
colnames(dm)=c("esd","model","w1","w2","w3")

g1<-ggplot(dm,aes(x=esd,y=w1,fill=model))+geom_bar(stat="identity")
g1=g1+theme_classic() # Clear background
g1=g1+scale_fill_manual(values=c("orange","purple")) # Set the colors
g1=g1+xlab("Error SD")+ylab("AIC weight") # Axis labels
g1=g1+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_mRNA_neutral.pdf",plot=g1,width=8,height=5) # Save the plot

g2<-ggplot(dm,aes(x=esd,y=w2,fill=model))+geom_bar(stat="identity")
g2=g2+theme_classic() # Clear background
g2=g2+scale_fill_manual(values=c("orange","purple")) # Set the colors
g2=g2+xlab("Error SD")+ylab("AIC weight") # Axis labels
g2=g2+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_translation_neutral.pdf",plot=g2,width=8,height=5) # Save the plot

g3<-ggplot(dm,aes(x=esd,y=w3,fill=model))+geom_bar(stat="identity")
g3=g3+theme_classic() # Clear background
g3=g3+scale_fill_manual(values=c("orange","purple")) # Set the colors
g3=g3+xlab("Error SD")+ylab("AIC weight") # Axis labels
g3=g3+theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=15)) # Adjust font and legend sizes
ggsave("aicw_protein_neutral.pdf",plot=g3,width=8,height=5) # Save the plot


