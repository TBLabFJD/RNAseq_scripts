#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=bowtieIndex   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gonzalo.n.moreno@gmail.com # Where to send mail        
##SBATCH --mem-per-cpu=10gb # Per processor memory
#SBATCH --cpus-per-task=80
#SBATCH -t 150:00:00     # Walltime
#SBATCH -o tophat2_%j.out # Name output file 
##SBATCH --error=
##SBATCH --file=
##SBATCH --initaldir=

# modules

module load tophat
module load bowtie2/2.3.4.3

# arguments

genome=$1
prefix=$2


# command line

bowtie2-build --threads 80 $genome $prefix 




