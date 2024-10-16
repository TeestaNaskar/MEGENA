#this script is for doing downstream data analysis of rat plecenta, all the scripts remain same for both vehicle and thc-cbd exposed group. Only the output directory will be chnaged as per the group.
setwd('/sc/arion/projects/MetaDope/Teesta/MEGENA/Rat.Placenta/')
library(MEGENA)
library(Matrix)
library(openxlsx)
library(matrixStats)
library(tidyverse)
library(dplyr)

#import data, normalize to counts per million
#IMPORT DATA
#import data, normalize to counts per million
#IMPORT DATA
data = read.csv("processed_counts.RatPlacenta/bothsex.rat.VST_counts.csv")

rownames(data)= toupper(data[,1])
data=data[,2:dim(data)[2]]

#IMPORT METADATA
meta=read.table("processed_counts.RatPlacenta/metadata.txt")
meta <- meta %>% rename (Group = Genotype)
meta <- meta %>% rename (Sex = Gender)

#Run MEGENA FOR EACH GROUP SEPERATELY (WRITE A LOOP)
#for (group in c("Cannabis", "Control")){
 # print(group)
  metasubset=meta[meta$Group== "THC_CBD",]
 animals_to_remove_from_metadata = setdiff(rownames(metasubset), colnames(data)) #identify animals to remove
   metasubset = metasubset[!(rownames(metasubset) %in% animals_to_remove_from_metadata),]

datExpr=data[,as.character(rownames(metasubset))]
  print(dim(datExpr))
setwd("CBD")
 #remove genes with stddev=0
    sd_rows <- apply(datExpr, 1, sd)
    datExpr = datExpr[which(sd_rows>0),]
    print(dim(datExpr))

set.seed(12345)

saveto="./output/"
dir.create(saveto)
sigmod=readLines("multiscale_significant.modules.txt")#significant modules#
sigmod=lapply(sigmod,function(x) strsplit(x,"\t")[[1]])
names(sigmod)=do.call(c,lapply(sigmod,function(x) x[1]))
sigmod=lapply(sigmod,function(x) x[-1])
sigmod=stack(sigmod)
write.table(sigmod,paste(saveto,"significant_module_2column_table.txt",sep=""),sep="\t",quote=F,row.names=F)


sigmod=readLines("multiscale_significant.modules.txt")#significant modules#
sigmod=lapply(sigmod,function(x) strsplit(x,"\t")[[1]])
names(sigmod)=do.call(c,lapply(sigmod,function(x) x[1]))
sigmod=lapply(sigmod,function(x) x[-1])
sigmod=stack(sigmod)
write.table(sigmod,paste(saveto,"significant_module_2column_table.txt",sep=""),sep="\t",quote=F,row.names=F)
