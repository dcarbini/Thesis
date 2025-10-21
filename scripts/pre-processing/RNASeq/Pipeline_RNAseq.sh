## Pipeline 

# In the server 
GSE='GSE153320' # modify it according to your dataset
cd /data/omics/public_data/endocrine_disruptors/
mkdir $GSE
mkdir $GSE/scripts $GSE/raw_data_GSM $GSE/raw_data $GSE/processed_data $GSE/raw_data_GSM_merged
mkdir $GSE/processed_data/bam $GSE/processed_data/counts $GSE/processed_data/trimmed_fastq $GSE/processed_data/unique_sorted

# Download the raw-data from ENA  

# Out of the server
file='ena-file-download-read_run-PRJNA642125-fastq_ftp-20250317-1505.sh' # modify it according to your dataset 
scp /mnt/c/Users/carbi/Downloads/$file hbp694@med-bart.rd.tuni.fi:/data/omics/public_data/endocrine_disruptors/$GSE/scripts/
#Put the password

# In the server 
tmux

cd /data/omics/public_data/endocrine_disruptors/$GSE/raw_data/
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/$file

#to enter in the conrol mode: ctrl+b 
#to exit: d 
#tmux ls
#tmux a 

# One time only:
#Download edirect or run the next 2 lines
cd
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo "export PATH=/home/hbp694/edirect:${PATH}" >> ${HOME}/.bashrc
export PATH=${HOME}/edirect:${PATH}

bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/SRR_GSM.sh $GSE
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/Merge_GSM.sh $GSE


# One time only:
source /data/miniforge3/bin/activate
conda init 
#log out and log in again

conda activate seqtools

cd /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/
fastqc /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/*
# Multiqc to combine multiple fastq reports in one overall report
multiqc . 

#Change the parameters in the script according the information retrieved from the fastqc reports before running it
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/trimming_alignment.sh $GSE

# Deactivate the old environment
conda deactivate
# Activate another environment to use R
conda activate renv

Rscript /data/omics/public_data/endocrine_disruptors/$GSE/scripts/count_matrix.R $GSE
