#this script is for megena analysis of Rat placenta samples specific for VEH group

setwd("/sc/arion/projects/MetaDope/Teesta/MEGENA/Rat.Placenta/")

library(MEGENA)
library(Matrix)
library(openxlsx)
library(matrixStats)
library(tidyverse)
options(warn=1)

#Megnea Input Parameters
n.cores=15
doPar=TRUE
method= "pearson" 
FDR.cutoff=2
module.pval=.05
hub.pval=.05
cor.perm=10 
hub.perm=10
min.size=10
max.size=NULL

#ANNOTATION TO BE DONE ON DOWNSTREAM
annot.table=NULL
id.col=1
symbol.col=2

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
  metasubset=meta[meta$Group== "VEH",]
 animals_to_remove_from_metadata = setdiff(rownames(metasubset), colnames(data)) #identify animals to remove
   metasubset = metasubset[!(rownames(metasubset) %in% animals_to_remove_from_metadata),]
  
datExpr=data[,as.character(rownames(metasubset))]
  print(dim(datExpr))
setwd("VEH")
 #remove genes with stddev=0
    sd_rows <- apply(datExpr, 1, sd)
    datExpr = datExpr[which(sd_rows>0),]
    print(dim(datExpr))

  #varianceofallgenes=rowVars(as.matrix(datasubset))
  #datasubset=datasubset[which(varianceofallgenes>0),]
  
  #if (group=="Cannabis"){setwd("Cannabis1")}
  #else if (group=="Control"){setwd("Control1")}
  
  #RUN ALL GENE CORRELATIONS
  ijw=calculate.correlation(datExpr=datExpr, doPerm=cor.perm, method=method, FDR.cutoff=FDR.cutoff, saveto=".")
  save(ijw,file="ijw.RData")
  print(dim(ijw))
  
  #SET UP CORES
  run.par=doPar & (getDoParWorkers()==1)
  if(run.par){
    cl=parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste("number of cores",getDoParWorkers(),"\n",sep = "" ))
  }
  #BUILD NETWORK
## Compute PFN
  el=calculate.PFN(ijw[,1:3], doPar=doPar,num.cores=n.cores)
  print("finished calculate PFN")
  
  #SAVE FILE
  write.csv(el,file="weights.for.cytoscape.csv", row.names=FALSE)
  saveRDS(el,file="PFN.RData")
  rm(ijw)
  
  #TAKING PFN AND TURNING IT INTO A GRAPH DATAFRAME, DATA STRCUTURE TO SPECIFY GENE EDGES/CONNECTIONS
  g=graph.data.frame(el,directed=FALSE)
  saveRDS(g,"Graph.RData")
  rm(el)
  
  #RUN MEGENA TO CLUSTER CORRELATIONS AND EDGES INTO MODULES
  megena.output=do.MEGENA(g,mod.pval=module.pval, hub.pval=hub.pval, remove.unsig=TRUE, min.size=min.size, 
                          max.size=vcount(g)/2,doPar=doPar, num.cores=n.cores,n.perm=hub.perm, save.output=T)
  print("Finished Megena")
  
  if(getDoParWorkers()>1){
    env=foreach:::.foreachGlobals
    rm(list = ls(name=env),pos=env)
  }
  
  #MAKE MODULE SUMMARY
  summary.output=MEGENA.ModuleSummary(megena.output, mod.pvalue=module.pval, hub.pvalue=hub.pval, 
                                      min.size=min.size, max.size=vcount(g)/2, annot.table=annot.table,id.col=id.col,
                                      symbol.col=symbol.col, output.sig=T)
  print("Finished Megena Module Summary")
  
  module.output=module_convert_to_table(megena.output,mod.pval=0.05,hub.pval=0.05, min.size=min.size, max.size=vcount(g)/2)
  save(summary.output, megena.output,module.output,g,file = "Megena.Results.RData")
  rm(list=c("summary.output","module.output"))
  setwd("../")
}
