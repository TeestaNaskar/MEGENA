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

library("msigdbi", lib.loc = "/sc/arion/projects/MetaDope/Teesta/R")
require(GOtest, lib.loc = "/sc/arion/projects/MetaDope/Teesta/R")
library(data.table)
sigmod = read.table("output/significant_module_2column_table.txt", header = T)


deg.total = read.xlsx("../processed_counts.RatPlacenta/bothsex.combined.TNdeseq_THC_CBD_vs_VEH.xlsx", sheet=1)
deg.total$Genes = toupper(deg.total$Genes)
query = deg.total$Genes

deg <- deg.total[deg.total$pvalue < 0.05 & abs(deg.total$log2FoldChange) >0,]
deg$Reg = ifelse (deg$log2FoldChange >0 , "UP", "DN")
deg_mod=GOtest(x=sigmod,deg[,c("Genes","Reg")],query.population = sigmod$values,background = "query",method="hypergeometric")
write.table(deg_mod,paste(saveto,"MEGENA_mod_overlap_DEGs.xlsx",sep=""),sep="\t",quote=F,row.names=F)
 msigdb.gsea -> msigdb.genesets -> match.arg
go=msigdb.gsea(sigmod,background = "annotation",method="hypergeometric",species = "rat")
write.table(go[go$P.adj<0.05,],paste(saveto,"MEGENA_mod_GO.xlsx",sep=""),sep="\t",quote=F,row.names=F)

tf=msigdb.gsea(sigmod,genesets = "c3.tft",background = "annotation",method = "hypergeometric",species = "rat")
write.table(tf[tf$P.adj<0.05,],paste(saveto,"MEGENA_mod_C3.TFT.xlsx",sep=""),sep="\t",quote=F,row.names=F)
