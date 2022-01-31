#!/bin/bash



module load python/2.7.15
module load perl
source ~/.Renviron



#output_path="/home/proyectos/bioinfo/fjd/reanalysis/rRiveiro/rRiveiro270121/CANCER31-09-10-2018/fastq"
#output_path="/home/proyectos/bioinfo/fjd/reanalysis/rRiveiro/rRiveiro270121/CANCER20-20-04-2017/fastq"
#pipeline="/home/proyectos/bioinfo/fjd/VariantCallingFJD/pipelineFJD19.py"
#output_path="/home/proyectos/bioinfo/fjd//reanalysis/mDelPozo/mdelpozo280121/fastq"
#output_path="/scratch/lodela/SureSelect-QXT1-20210126"

#/home/proyectos/bioinfo/software/bs download project --name "CANCER31-09-10-2018" -o $output_path --extension=fastq.gz

#/home/proyectos/bioinfo/software/bs download project --name "CANCER20-20-04-2017" -o $output_path --extension=fastq.gz

#/home/proyectos/bioinfo/software/bs download project --name "SureSelect-QXT1-20210126" -o $output_path --extension=fastq.gz
#output_path="/scratch/lodela/20210209_Carrera_SureselectQXT_2"
#/home/proyectos/bioinfo/software/bs download project --name "20210209_Carrera SureselectQXT_2" -o $output_path --extension=fastq.gz

#/home/proyectos/bioinfo/software/bs download run --name  "SureSelect-QXT1-20210126" -o $output_path
#batch --account=bioinfo_serv --partition=bioinfo  /home/proyectos/bioinfo/lodela/BioinfoUnit/slurm_scripts/bcltofastq.sh $manifest ./SureSelect-QXT1-20210126/ SureSelect-QXT1-20210126_fastq

#output_path="/scratch/lodela/20210209_Carrera_SureselectQXT_2"
#/home/proyectos/bioinfo/software/bs download run --name  "20210209_Carrera SureselectQXT_2" -o $output_path
#sbatch --account=bioinfo_serv --partition=bioinfo  /home/proyectos/bioinfo/lodela/BioinfoUnit/slurm_scripts/bcltofastq.sh $manifest ./SureSelect-QXT1-20210126/ SureSelect-QXT1-20210126_fastq

output_path="/scratch/lodela/SureSelect-QXT3-20210301"
#/home/proyectos/bioinfo/software/bs download run --name  "SureSelect-QXT3-20210301_ERROR nºADNs!!!!!" -o $output_path
manifest="./SureSelect-QXT3-20210301/sample_sheet_QXT3_ADNcorrecto_29-03-2021.csv"

sbatch --account=bioinfo_serv --partition=bioinfo  /home/proyectos/bioinfo/lodela/BioinfoUnit/slurm_scripts/bcltofastq.sh $manifest $output_path ${output_path}ADNcorrecto_fastq
echo sbatch --account=bioinfo_serv --partition=bioinfo  /home/proyectos/bioinfo/lodela/BioinfoUnit/slurm_scripts/bcltofastq.sh $manifest $output_path ${output_path}ADNcorrecto_fastq


#/home/proyectos/bioinfo/software/bs-cp //mRodriguez/Run/20200120Panic Linfomas_NextSeq500_200420/




#output_path="/home/proyectos/bioinfo/NOBACKUP/lodela/EPIC"
#/home/proyectos/bioinfo/software/bs download project --name "EPIC4-2021-02-01" -o $output_path --extension=fastq.gz




#output_path="/home/proyectos/bioinfo/NOBACKUP/reanalysis2/rRiveiro/rRiveiro50321/fastq"
#/home/proyectos/bioinfo/software/bs download project --name "CANCER34-20-02-19" -o $output_path --extension=fastq.gz
