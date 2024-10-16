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
  #metasubset=meta[meta$Group== "THC_CBD",]
 animals_to_remove_from_metadata = setdiff(rownames(meta), colnames(data)) #identify animals to remove
   metasubset = meta[!(rownames(meta) %in% animals_to_remove_from_metadata),]

datExpr=data[,as.character(rownames(metasubset))]
  print(dim(datExpr))
setwd("allsamples")
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
#go=msigdb.gsea(sigmod,background = "annotation",method="hypergeometric",species = "rat")
#write.table(go[go$P.adj<0.05,],paste(saveto,"MEGENA_mod_GO.xlsx",sep=""),sep="\t",quote=F,row.names=F)

#tf=msigdb.gsea(sigmod,genesets = "c3.tft",background = "annotation",method = "hypergeometric",species = "rat")
#write.table(tf[tf$P.adj<0.05,],paste(saveto,"MEGENA_mod_C3.TFT.xlsx",sep=""),sep="\t",quote=F,row.names=F)


library(DGCA, lib = "/sc/arion/projects/MetaDope/Teesta/R")

design_mat=model.matrix(~0+meta$Group)
head(design_mat)

 colnames(design_mat)=c("THC_CBD","VEH")
moduleDC_res=moduleDC(inputMat = datExpr,design = design_mat,
                      compare = c("VEH","THC_CBD"),genes=sigmod$values,
                      labels=sigmod$ind,nPerm=50, gene_avg_signif = 1,
                      number_DC_genes = 100)
write.table(moduleDC_res,paste0(saveto,"ds_MDC_nperm25.xlsx"),sep="\t",quote=F,row.names = F)


load("../MEGENA.Results.RData")
#visualization of pnet object
pnet.obj <- plot_module(output.summary = summary.output,PFN = g, layout = "kamada.kawai",
                        label.hubs.only = FALSE,
                        gene.set = NULL,
                        output.plot = TRUE,
                        out.dir = "module_plot",
                        col.names = c("red", "orange", "yellow", "green", "blue", "purple"),
                        label.scaleFactor = 40,
                        hubLabel.col = "black",
                        hubLabel.sizeProp = 1,
                        show.topn.hubs = Inf,
                        show.legend = TRUE)
saveRDS(pnet.obj, file="pnet_obj.RData")

#X11();
print(pnet.obj[[1]])

module.table <- summary.output$module.table
colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".

hierarchy.obj <- plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
                                       arrow.size = 0.03,node.label.color = "blue")
#X11();
print(hierarchy.obj[[1]])

saveRDS(hierarchy.obj, file="hierarchy.RData")


topo=MEGENA.output$module.output$module.relation
head(topo)
topo[,1]=paste0("c1_",topo[,1])
topo[,2]=paste0("c1_",topo[,2])
colnames(topo)=c("mod.parent","module.id")
topo=as.data.frame(topo)


get.topEnrich=function(enrichtable,subject,byfactor){
  s=split.data.frame(enrichtable,as.factor(enrichtable[,match(subject,colnames(enrichtable))]))
  ss=lapply(s,function(x,byfactor) {
    a=as.numeric(as.character(x[,match(byfactor,colnames(x))]))
  index=ifelse(min(a)<0.05,which.min(a),NA);
  sss=do.call(rbind,ss)
  sss=sss[!is.na(sss[,1]),]
  return(sss)
  })}


deg_mod1=get.topEnrich(deg_mod,subject = "Input",byfactor = "P.adj")
deg_mod1$log10FDR=ifelse(deg_mod1$Category=="UP",-log10(deg_mod1$P.adj),log10(deg_mod1$P.adj))
topo$DEG=deg_mod1$log10FDR[match(topo[,2],deg_mod1$Input)]

ct1=get.topEnrich(ct,subject = "Input",byfactor = "P.adj")
topo$celltype=ct1$Category[match(topo$module.id,ct1$Input)]

md=moduleDC_res[moduleDC_res$pVal<0.05,]
topo$MDC=md$MeDC[match(topo$module.id,md$Module)]

size=as.data.frame.vector(table(sigmod$ind))
topo$size=size[match(topo$module.id,rownames(size)),1]

go1=go[go$P.adj<0.05,]
go1=get.topEnrich(go1,subject = "Input",byfactor = "P.adj")
topo$topGO=go1$MSigDB[match(topo$module.id,go1$Input)]

tf1=tf[tf$P.adj<0.05,]
tf1=get.topEnrich(tf1,subject = "Input",byfactor = "P.adj")
topo$topTF=tf1$MSigDB[match(topo$module.id,tf1$Input)]

topo[is.na(topo)]=0
topo$Ranking=2*nrow(topo)-rank(abs(topo$DEG),ties.method = "average")-rank(abs(topo$MDC),ties.method = "average")
topo=topo[order(topo$Ranking,decreasing = F),]
topo[topo==0]=NA
write.table(topo,paste0(saveto,"module.summary.txt"),sep="\t",quote=F,row.names = F)
