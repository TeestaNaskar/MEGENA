#this script is for the downstream analysis enriching with male rat DEGs with created MEGENA co-expression network for all samples. 
#So, first MEGENA ran with all samples to create network modules and then enriched particularly with DEGs only for males and run GO on that 

setwd('/sc/arion/projects/MetaDope/Teesta/MEGENA/Rat.Placenta/allsamples/males')
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
data = read.csv("../../processed_counts.RatPlacenta/bothsex.rat.VST_counts.csv")

rownames(data)= toupper(data[,1])
data=data[,2:dim(data)[2]]

#IMPORT METADATA

meta = read.table("../../processed_counts.RatPlacenta/metadata.txt")

meta <- meta %>% rename (Group = Genotype)
meta <- meta %>% rename (Sex = Gender)

#Run MEGENA FOR EACH GROUP SEPERATELY (WRITE A LOOP)
#for (group in c("Cannabis", "Control")){
 # print(group)
  #metasubset=meta[meta$Group== "THC_CBD",]
#remove X from the rownames
#rownames(meta) = gsub('X', "", rownames(meta))

animals_to_remove_from_metadata = setdiff(rownames(meta), colnames(data)) #identify animals to remove
   metasubset = meta[!(rownames(meta) %in% animals_to_remove_from_metadata),]

datExpr=data[,as.character(rownames(metasubset))]
  print(dim(datExpr))

saveto="./output/"
dir.create(saveto)

#remove genes with stddev=0
    sd_rows <- apply(datExpr, 1, sd)
    datExpr = datExpr[which(sd_rows>0),]
    print(dim(datExpr))

set.seed(12345)


sigmod=readLines("../combined_sex/multiscale_significant.modules.txt")#significant modules#
sigmod=lapply(sigmod,function(x) strsplit(x,"\t")[[1]])
names(sigmod)=do.call(c,lapply(sigmod,function(x) x[1]))
sigmod=lapply(sigmod,function(x) x[-1])
sigmod=stack(sigmod)

# Convert Gene Symbols to Title Case by following the rat/mouse nomenclature
library(stringr)
sigmod$values <- str_to_title(sigmod$values)  # "Apoc2", "Fgg", etc.

write.table(
  sigmod,
  file = paste0(saveto, "significant_module_2column_table.txt"),
  sep = "\t",          # Use proper tab separator
  quote = FALSE,       # Prevent quotes around values
  row.names = FALSE    # Omit row names
)

library("msigdbi", lib.loc = "/sc/arion/projects/MetaDope/Teesta/R")
require(GOtest, lib.loc = "/sc/arion/projects/MetaDope/Teesta/R")
library(data.table)
sigmod <- read.table(
  "output/significant_module_2column_table.txt",
  header = TRUE,       # Read the first row as column headers
  sep = "\t",          # Specify tab as the separator
  stringsAsFactors = FALSE  # Ensure character columns aren't converted to factors
)

#loading the male DEG file for enriching

deg = read.csv("../../processed_counts.RatPlacenta/DEG/teesta.deseq_THC_CBD_Male_vs_VEH_Male.csv")

#deg$Genes = toupper(deg$Genes)
deg$Genes <- str_to_title(deg$Genes) #to confrim every degs are in same case
query = deg$Genes

deg <- deg[deg$pvalue < 0.05 & abs(deg$log2FoldChange) >0,]
deg$Reg = ifelse (deg$log2FoldChange >0 , "UP", "DN")

deg_mod=GOtest(x=sigmod,deg[,c("Genes","Reg")],query.population = sigmod$values,background = "query",method="hypergeometric")

write.xlsx(deg_mod,"output/MEGENA_mod_overlap_DEGs.xlsx")

go = msigdb.gsea(sigmod,background = "annotation",method="hypergeometric",species = "mouse")
write.xlsx(go[go$P.adj<0.05,],"output/MEGENA_mod_GO.xlsx")

tf=msigdb.gsea(sigmod,genesets = "c3.tft",background = "annotation",method = "hypergeometric",species = "mouse")

write.xlsx(tf[tf$P.adj<0.05,],"output/MEGENA_mod_C3.TFT.xlsx")


load("../combined_sex/Megena.Results.RData")

topo=megena.output$module.output$module.relation
head(topo)
topo[,1]=paste0("c1_",topo[,1])
topo[,2]=paste0("c1_",topo[,2])
colnames(topo)=c("mod.parent","module.id")
topo=as.data.frame(topo)

get.topEnrich <- function(enrichtable, subject, byfactor) {
  # Split the table based on the "subject" column
  s <- split(enrichtable, as.factor(enrichtable[, match(subject, colnames(enrichtable))]))
  
  # Apply function to each split table
  ss <- lapply(s, function(x) {
    # Extract the column specified by "byfactor"
    a <- as.numeric(as.character(x[, match(byfactor, colnames(x))]))
    
    # Find the index of the minimum value if it's below 0.05
    index <- ifelse(min(a, na.rm = TRUE) < 0.05, which.min(a), NA)
    
    # Return the corresponding row if an index was found
    if (!is.na(index)) {
      return(x[index, , drop = FALSE])  # Ensure row is returned as a data frame
    } else {
      return(NULL)
    }
  })
  
  # Combine results into a single data frame
  sss <- do.call(rbind, ss)
  
  # Remove rows with NA results
  sss <- sss[!is.na(sss[, 1]), ]
  
  return(sss)
}


deg_mod1=get.topEnrich(deg_mod,subject = "Input",byfactor = "P.adj")
deg_mod1$log10FDR=ifelse(deg_mod1$Category=="UP",-log10(deg_mod1$P.adj),log10(deg_mod1$P.adj))
topo$DEG=deg_mod1$log10FDR[match(topo[,2],deg_mod1$Input)]

#ct1=get.topEnrich(ct,subject = "Input",byfactor = "P.adj")
#topo$celltype=ct1$Category[match(topo$module.id,ct1$Input)]

#md=moduleDC_res[moduleDC_res$pVal<0.05,]
#topo$MDC=md$MeDC[match(topo$module.id,md$Module)]

size=as.data.frame.vector(table(sigmod$ind))
topo$size=size[match(topo$module.id,rownames(size)),1]

go1=go[go$P.adj<0.05,]
go1=get.topEnrich(go1,subject = "Input",byfactor = "P.adj")
topo$topGO=go1$MSigDB[match(topo$module.id,go1$Input)]

tf1=tf[tf$P.adj<0.05,]
tf1=get.topEnrich(tf1,subject = "Input",byfactor = "P.adj")
topo$topTF=tf1$MSigDB[match(topo$module.id,tf1$Input)]

topo[is.na(topo)]=0
#topo$Ranking=2*nrow(topo)-rank(abs(topo$DEG),ties.method = "average")-rank(abs(topo$MDC),ties.method = "average")
#topo=topo[order(topo$Ranking,decreasing = F),]
topo[topo==0]=NA
write.table(topo,paste0(saveto,"module.summary.txt"),sep="\\t",quote=F,row.names = F)
