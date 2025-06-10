#! /bin/bash

#SBATCH --job-name=lmt_preprocess_align
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL 
#SBATCH --mail-user=cesar.martinez@isglobal.org 
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100gb
#SBATCH --output=/PROJECTES/MALARIA_EPIGENETICS/Job_Logs/cmartinez_%j.log

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the bashrc and conda environment 
module load lang/Anaconda3/2023.03-1
source ~/.bashrc
source activate env37


# And finally run the job
#/home/isglobal.lan/lmichel/.conda/envs/python36/bin/python3.6 /PROJECTES/MALARIA_EPIGENETICS/Scripts/Pipelines/chip_seq_pipeline.py -pf /PROJECTES/MALARIA_EPIGENETICS/Projects/Alba/Reanalysis_C2/chipseq_pipeline_parameters_C2.txt
#/home/isglobal.lan/lmichel/.conda/envs/lucas_main/bin/python3.8 /PROJECTES/MALARIA_EPIGENETICS/Scripts/Pipelines/chip_seq_pipeline.py -pf /PROJECTES/MALARIA_EPIGENETICS/Projects/PhD_Project/chipseq_pipeline_parameters_12B.txt
#/home/isglobal.lan/cmartinez/.conda/envs/cesar_main/bin/python3.11 /PROJECTES/MALARIA_EPIGENETICS/Scripts/Pipelines/chip_seq_pipeline.py -pf /PROJECTES/MALARIA_EPIGENETICS/ChIP_Seqs/ETF_33-34/chipseq_pipeline_parameters.txt
python3 /PROJECTES/MALARIA_EPIGENETICS/Scripts/Pipelines/Pipeline_Genewise_Coverage/genewise_coverage_pipeline.py -pf /PROJECTES/MALARIA_EPIGENETICS/Scripts/Pipelines/Pipeline_Genewise_Coverage/genewise_coverage_pipeline_parameters.txt
