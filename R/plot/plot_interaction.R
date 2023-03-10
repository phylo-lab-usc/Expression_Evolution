#setwd("/Users/daohanji/Desktop/Expression_Evolution")
setwd("/Users/rexjiang/Desktop/Lab/multivariate_trait/expression/rerun/out_interaction")

library(ggplot2)

pt=1 # Part 1 = only one gene is under direct selection, part 2 = both genes are inder direct selection

# Input file, can be out_interaction_2g_pt1_end.txt or out_interaction_2g_pt2_end.txt
fn=paste("out_interaction_2g_pt",pt,"_end.txt",sep="")
d<-read.table(fn,sep="\t")
# Set column names
colnames(d)=c("c1","c2",
	"var_transcription1",
	"var_translation1",
	"var_transcription2",
	"var_translation2",
	"var_mRNA1",
	"var_protein1",
	"var_mRNA2",
	"var_protein2",
	"cor_gv1",
	"cor_gv2",
	"cor_phe1",
	"cor_phe2",
	"var_ratio1",
	"var_ratio2")

# Make heatmap
h1<-ggplot(d,aes(x=c1,y=c2,fill=cor_gv1))+geom_tile()
h1=h1+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL) # Set the color scheme and remove legend text (to be described in the figure legend instead)
h1=h1+theme_classic() # Clear background
h1=h1+theme(axis.text=element_text(size=15),legend.text=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust text size
h1=h1+xlab("")+ylab("") # Remove axis labels, to be custom added later

h2<-ggplot(d,aes(x=c1,y=c2,fill=cor_gv2))+geom_tile()
h2=h2+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL)
h2=h2+theme_classic() # Clear background
h2=h2+theme(axis.text=element_text(size=15),legend.text=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust text size
h2=h2+xlab("")+ylab("") # Remove axis labels, to be custom added later

h3<-ggplot(d,aes(x=c1,y=c2,fill=cor_phe1))+geom_tile()
h3=h3+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL) # Set the color scheme and remove legend text (to be described in the figure legend instead)
h3=h3+theme_classic() # Clear background
h3=h3+theme(axis.text=element_text(size=15),legend.text=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust text size
h3=h3+xlab("")+ylab("") # Remove axis labels, to be custom added later

h4<-ggplot(d,aes(x=c1,y=c2,fill=cor_phe2))+geom_tile()
h4=h4+scale_fill_gradient2(low="blue4",mid="white",high="gold",limits=c(-1,1),name=NULL) # Set the color scheme and remove legend text (to be described in the figure legend instead)
h4=h4+theme_classic() # Clear background
h4=h4+theme(axis.text=element_text(size=15),legend.text=element_text(size=15),legend.key.size=unit(0.7,'cm')) # Adjust text size
h4=h4+xlab("")+ylab("") # Remove axis labels, to be custom added later

# Save the plots
fn1=paste("interaction_2g_pt",pt,"_cor1_g1.pdf",sep="")
fn2=paste("interaction_2g_pt",pt,"_cor1_g2.pdf",sep="")
fn3=paste("interaction_2g_pt",pt,"_cor2_g1.pdf",sep="")
fn4=paste("interaction_2g_pt",pt,"_cor2_g2.pdf",sep="")
ggsave(fn1,plot=h1,width=6.25,height=5)
ggsave(fn2,plot=h2,width=6.25,height=5)
ggsave(fn3,plot=h3,width=6.25,height=5)
ggsave(fn4,plot=h4,width=6.25,height=5)


