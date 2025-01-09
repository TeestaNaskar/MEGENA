setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Rat_placenta/RNAseq/MEGENA/Rat.Placenta/THC_CBD")
library(data.table)
library(stringr)
NodeSum=fread("multiscale_nodeSummary.txt")
head(NodeSum)
NodeSum[1:5,1:5]

# Providing the number modules the gene is a key driver in
KDlength=c()
for (KD in 1:nrow(NodeSum)){
  mods=NodeSum$KD.membership[KD]
  if (nchar(mods)==0){
    KDlength=c(KDlength,0)
    next
  }
  KDlength=c(KDlength,length(str_split(mods,",")[[1]]))    
}

NodeSum$nummodsKD=KDlength
NodeSum=NodeSum[order(-nummodsKD),]

trimmednodsum=NodeSum[,c("id","module.membership","KD.membership","nummodsKD")]
write.csv(trimmednodsum,file="keydrivers_THC.Rats.csv")
