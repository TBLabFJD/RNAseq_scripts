#!/bin/bash

index="/home/proyectos/bioinfo/references/bowtie_index_hg38/hg38"
gtf="/home/proyectos/bioinfo/fjd/BioinfoUnit/gffs/hsa_maturemiRNA.gff3"

output_dir="/home/proyectos/bioinfo/NOBACKUP/BioinfoUnit/onco"

script_dir="${output_dir}/scripts"
fastq_dir="${output_dir}/fastq"
log_dir="${output_dir}/logfiles"

mkdir ${output_dir}/trimmed_fastq/
mkdir ${output_dir}/mapped/
mkdir ${output_dir}/counts/
mkdir ${output_dir}/qc/
mkdir ${log_dir}






###################
# Quality control #
###################

for file in ${fastq_dir}/*.gz
do 
sample=$(basename "$file")
sbatch --account=bioinfo_serv --partition=fastbioinfo \
-e ${log_dir}/${sample}_qc.err -o ${log_dir}/${sample}_qc.out \
${script_dir}/fastqc_slurm.sh ${file} ${output_dir}/qc/ 
done






############
# trimming #
############
for file in ${fastq_dir}/*R1.f*q.gz
do

# NEXTflex
sbatch --account=bioinfo_serv --partition=bioinfo \
-e ${log_dir}/${sample}_cutadapt.err -o ${log_dir}/${sample}_cutadapt.out \
${script_dir}/cutadapt_NEXTflex_slurm.sh ${file} ${output_dir}/trimmed_fastq/

# # TruSeq
# sbatch --account=bioinfo_serv --partition=bioinfo \
# -e ${log_dir}/${sample}_cutadapt.err -o ${log_dir}/${sample}_cutadapt.out \
# ${script_dir}/cutadapt_TruSeq_slurm.sh ${file} ${output_dir}/trimmed_fastq/

done






###########
# Mapping #
###########

for file in ${output_dir}/trimmed_fastq/*_l17_25.fq
do
sample=$(basename "$file" _l17_25.fq)
echo ${sample}


# miRNA R1 y R2
job=`sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=16 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_bowtie.err -o ${log_dir}/${sample}_bowtie.out \
${script_dir}/bowtie_slurm.sh ${index}  ${R1} ${R2} ${output_dir}/mapped/`


jobid=$(echo $job | cut -f4 -d" ")
bam=${output_dir}/mapped/${sample}.bam
output_count=${output_dir}/counts/${sample}.counts.qualityUniqueStranded.txt

#sbatch --account=bioinfo_serv --partition=bioinfo -e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out --dependency=afterok:$jobid ${script_dir}/htseqMiRNA_slurm.sh $bam $gtf ${output}
sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=20 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out --dependency=afterok:$jobid \
${script_dir}/htseqMiRNA_slurm.sh $bam $gtf ${output_count}

done






#########
# Count #
#########

for bam in ${output_dir}/mapped/*
do
sample=$(basename "$bam" accepted_hits.bam)
echo ${sample}

sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=20 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out \
${script_dir}/htseqMiRNA_slurm.sh $bam $gtf ${output_count}

done





