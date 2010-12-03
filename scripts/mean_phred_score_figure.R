#!/usr/bin/env Rscript
# Takes an error profile as an argument, then loads it and produces a
# boxplot figure showing the phred scores by position and reference base
#
#Call with `Rscript boxplot_phred_score_figure.R [input filename] [output filename]
Args <- commandArgs(TRUE) #grab command args
if(length(Args) != 2){
  print("Usage: Rscript mean_phred_score_figure.R input.txt output.[png,pdf,tiff,tex,ps,...]")
  quit()
}
library(reshape)
library(ggplot2)
tmp.dat=read.table(Args[0])
tmp.df=data.frame(tmp.dat)
names(tmp.df)<-c("Position", "Phred","A->A","A->C","A->G","A->T","A->N","C->A","C->C","C->G","C->T","C->N","G->A","G->C","G->G","G->T","G->N","T->A","T->C","T->G","T->T","T->N")
tmp.melted<-melt(tmp.df,id.var=c("Position","Phred"))
names(tmp.melted)<-c("Position","Phred","Substitution","Frequency")
tmp.melted$Reference<-substr(tmp.melted$Substitution,1,1)
qplot(Position,Phred,data=tmp.melted,weight=Frequency,color=Substitution,geom="smooth",facets=.~ Reference)
ggsave(filename=Args[2])

