#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=fastqc   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gonzalo.n.moreno@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=2
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_fastqc.out # Name output file 
#SBATCH --error=%j_fastqc.err
##SBATCH --file=
##SBATCH --initaldir=

# modules
module load fastqc/0.11.7

# arguments

fastq=$1
output_dir=$2

# command line

fastqc $fastq -t 2 -o $output_dir
