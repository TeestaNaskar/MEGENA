#BSUB -L /bin/bash
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 20
#BSUB -R span[ptile=14]
#BSUB -R rusage[mem=30000]
#BSUB -W 30:00
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I

module load R/4.0.3
module load openssl
module load udunits


R CMD BATCH downstream.analysis.allsamples.R 
##the R script will be changed as per the name of the R script prepared for doing downstream analysis. it might be complete downstream R script or 
