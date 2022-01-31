#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=htseqcount   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=gonzalo.n.moreno@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=40
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_htseq.out # Name output file 
#SBATCH --error=%j_htseq.err
##SBATCH --file=
##SBATCH --initaldir=

# modules
module load miniconda-bio/3.7

# arguments

bam=$1
gtf=$2
output=$3

echo $bam
echo $gtf
echo $output


# command line

#python $htseq_count $bam $gtf -f bam -r pos -t exon -i gene_id --additional-attr=gene_name --stranded=reverse > ${output}


htseq-count $bam $gtf --stranded=no -a 0 -f bam -r pos -t exon -i gene_id --additional-attr=gene_name --max-reads-in-buffer=60000000 > ${output}






