#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=tophat2   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gonzalo.n.moreno@gmail.com # Where to send mail        
##SBATCH --mem-per-cpu=7gb # Per processor memory
#SBATCH --cpus-per-task=16
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_tophat2.out # Name output file 
#SBATCH --error=%j_tophat2.err
##SBATCH --file=
##SBATCH --initaldir=

# modules

module load samtools/1.9
module load tophat
module load python/2.7.15

# arguments

output=$1
bowtie_index=$2
R1=$3

echo $bowtie_index
echo $R1
echo $output



# command line

mkdir $output

tophat2 -p 15 -o $output $bowtie_index $R1







