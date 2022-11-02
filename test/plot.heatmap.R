setwd("/Users/daohanji/Desktop/Expression_Evolution")

library(ggplot2)

# Input file, can be out_interaction_2g_pt1_end.txt or out_interaction_2g_pt2_end.txt
d<-read.table("out_interaction_2g_pt2_end.txt",sep="\t")
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
h1=h1+scale_fill_gradient(low="blue",high="red",limits=c(-1,1))+theme_classic()

h2<-ggplot(d,aes(x=c1,y=c2,fill=cor_gv2))+geom_tile()
h2=h2+scale_fill_gradient(low="blue",high="red",limits=c(-1,1))+theme_classic()

h3<-ggplot(d,aes(x=c1,y=c2,fill=var_ratio1))+geom_tile()
h3=h3+scale_fill_gradient(low="purple",high="yellow",limits=c(0,3))+theme_classic()

h4<-ggplot(d,aes(x=c1,y=c2,fill=var_ratio2))+geom_tile()
h4=h4+scale_fill_gradient(low="purple",high="yellow",limits=c(0,3))+theme_classic()