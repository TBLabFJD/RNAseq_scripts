#!/bin/bash

# Human
index="/home/proyectos/bioinfo/references/bowtie_index_hg38/hg38" 
gtf="/home/proyectos/bioinfo/fjd/BioinfoUnit/gffs/gencode.v38.annotation.gff3"

# Mouse
index="/home/proyectos/bioinfo/references/bowtie_index_mm10/mm10"
gtf="/home/proyectos/bioinfo/fjd/BioinfoUnit/gffs/gencode.v38.annotation.gff3"


output_dir="/home/proyectos/bioinfo/NOBACKUP/BioinfoUnit/onco"

script_dir="${output_dir}/scripts"
fastq_dir="${output_dir}/fastq"
log_dir="${output_dir}/logfiles"

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
${script_dir}/fastqc_slurm.sh $file ${output_dir}/qc/ 
done






#####################
# Mapping and count #
#####################

for file in ${fastq_dir}/*R1.f*q.gz
do
sample=$(basename "$file" | sed 's/_R1.*//')

echo ${sample}


# RNAseq R1 y R2
R1=${output_dir}/fastq/${sample}_R1.fq.gz
R2=${output_dir}/fastq/${sample}_R2.fq.gz
job=`sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=16 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_tophat.err -o ${log_dir}/${sample}_tophat.out \
${script_dir}/tophat_slurm.sh ${output_dir}/mapped/${sample} ${index} ${R1} ${R2}`

# RNAseq R1
#job=`sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=16 --mail-user=gonzalo.n.moreno@gmail.com \
#-e ${log_dir}/${sample}_tophat.err -o ${log_dir}/${sample}_tophat.out \
#${script_dir}/tophatR1_slurm.sh ${output_dir}/mapped/${sample} ${index} ${file}`



jobid=$(echo $job | cut -f4 -d" ")
bam=${output_dir}/mapped/${sample}/accepted_hits.bam
output_count=${output_dir}/counts/${sample}.counts.qualityUniqueStranded.txt

#sbatch --account=bioinfo_serv --partition=bioinfo -e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out --dependency=afterok:$jobid ${script_dir}/htseqMiRNA_slurm.sh $bam $gtf ${output}
sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=20 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out --dependency=afterok:$jobid \
${script_dir}/htseq_slurm.sh ${bam} ${gtf} ${output_count}

done






#########
# Count #
#########

for bam_dir in ${output_dir}/mapped/*
do

sample=$(basename "${bam_dir}" accepted_hits.bam)
bam="${bam_dir}/accepted_hits.bam"
echo ${sample}
sbatch --account=bioinfo_serv --partition=bioinfo --cpus-per-task=20 --mail-user=gonzalo.n.moreno@gmail.com \
-e ${log_dir}/${sample}_count.err -o ${log_dir}/${sample}_count.out \
${script_dir}/htseq_slurm.sh ${bam} ${gtf} ${output_count}

done





