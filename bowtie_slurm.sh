#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=fastbioinfo
#SBATCH --job-name=bowtie2   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gonzalo.n.moreno@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=2
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_bowtie2.out # Name output file 
#SBATCH --error=%j_bowtie2.err
##SBATCH --file=
##SBATCH --initaldir=

# modules

module load samtools/1.9
module load bowtie2/2.3.4.3
module load python/2.7.15

# arguments

bowtie_index=$1
R1=$2
R2=$3
output=$4

echo $bowtie_index
echo $R1
echo $R2
echo $output


# command line: mapping plus samtobam

mkdir $output

sample=$(basename "$R1" .gz)
sample=$(basename "$sample" .fq)
sample=$(basename "$sample" .fastq)
sample=$(basename "$sample" _R1)


sam=${output}/${sample}.sam
bam=${output}/${sample}.bam

echo $sample
echo $bam
echo $sam


#bowtie2 -p 2 --very-sensitive-local -x $bowtie_index $R1  | samtools view -bS  | samtools sort -n --threads 2 -O bam -o $bam
bowtie2 -p 2  -x $bowtie_index $R1 $R2  | samtools view -bS  | samtools sort -n --threads 2 -O bam -o $bam
