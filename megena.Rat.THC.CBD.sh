#BSUB -L /bin/bash 
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 15
#BSUB -R span[ptile=14]
#BSUB -R rusage[mem=45000]
#BSUB -W 100:00
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I

module load R/4.0.3
module load openssl
module load udunits


R CMD BATCH megena.totalgroup.Rat.Placenta.R #this R script will get changed for the R script we want to run for eg. VEH or THC-CBD or for total group
