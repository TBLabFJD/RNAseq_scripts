#!/bin/sh

##### MiRNA mapping and quantification


#Pipeline:
#	-Mapping -> bowtie2
#	-Quantification --> htseq-count


# variables:


sampleDir="/home/proyectos/bioinfo/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_28Dec/trimmed_data_nextflex_l17_25"
mappedDir="/home/proyectos/bioinfo/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_28Dec/mapped_data_nextflex_l17_25"
index="/home/proyectos/bioinfo/references/bowtie_index_hg38/hg38" 
countsDir="/home/proyectos/bioinfo/fjd/BioinfoUnit/NataliaMarta_Vito/Marta_28Dec/counts_data_wAMB_stranded_nextflex_l17_25"

gtf="/home/proyectos/bioinfo/lodela/BioinfoUnit/NataliaVito/hsa_maturemiRNA.gff3"
#genome="/home/proyectos/bioinfo/references/hg38/hg38.fa.gz"

mkdir $mappedDir
mkdir $countsDir



# slurm:

bowtie="/home/proyectos/bioinfo/fjd/BioinfoUnit/slurm_scripts/bowtie_slurm.sh" 
htseq_count="/home/proyectos/bioinfo/fjd/BioinfoUnit/slurm_scripts/htseqMiRNA_slurm.sh" 

#bowtie2_build="/home/proyectos/bioinfo/lodela/BioinfoUnit/slurm_scripts/bowtie2build_slurm.sh" 


mail="--mail-user=gonzalo.n.moreno@gmail.com"



for file in $sampleDir/*_l17_25.fq; do


	echo $file
	sample=$(basename "$file" .fq)
	echo $sample

	# MAPPING
	job=`sbatch -e ${sample}_bowtie.err -o ${sample}_bowtie.out --account=bioinfo_serv --partition=bioinfo $mail $bowtie ${index} $file $mappedDir`
	jobid=$(echo $job | cut -f4 -d" ")
	bam=${mappedDir}/${sample}.bam
	echo $bam
	output=${countsDir}/${sample}.counts.qualityUniqueStranded.txt

	echo $output

	# READ EXTRACTION
	sbatch --account=bioinfo_serv --partition=bioinfo -e ${sample}_htseq.err -o ${sample}_htseq.out  --dependency=afterok:$jobid $mail $htseq_count $bam $gtf $output


			
done


