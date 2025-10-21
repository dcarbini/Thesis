# RNA Seq pipeline
Write a pipeline to pre-process RNA-seq data.

## 0 - Prepare the data structure and upload the metadata

``` 
# In the server 
GSE='GSE153320' # modify it according to your dataset
cd /data/omics/public_data/endocrine_disruptors/
mkdir $GSE
mkdir $GSE/scripts $GSE/raw_data_GSM $GSE/raw_data $GSE/processed_data $GSE/raw_data_GSM_merged $GSE/results $GSE/metadata
mkdir $GSE/processed_data/bam $GSE/processed_data/counts $GSE/processed_data/trimmed_fastq $GSE/processed_data/unique_sorted
```

```
# Out of the server
cd /mnt/c/Users/carbi/Desktop/UNIVERSITA/MAGISTRALE/BOLOGNA/Erasmus/Tampere/Thesis/esperanto/CURATION/$GSE/
scp $GSE.xlsx hbp694@med-bart.rd.tuni.fi:/data/omics/public_data/endocrine_disruptors/$GSE/metadata/
#Put the password
scp last_round_updated_$GSE.xlsx hbp694@med-bart.rd.tuni.fi:/data/omics/public_data/endocrine_disruptors/$GSE/metadata/updated_$GSE.xlsx
#or
scp updated_$GSE.xlsx hbp694@med-bart.rd.tuni.fi:/data/omics/public_data/endocrine_disruptors/$GSE/metadata/
#Put the password
```

## 1 - Download the raw-data
The raw data will be downloaded from ENA (https://www.ebi.ac.uk/ena/browser/home).
* Search the GSE in the box 
* Select the study and show the Read Files 
* Download All the FASTQ files 
	** NOTE: 
	1. Check if the dataset is TempO-Seq. If yes, skip all this part and start after point 4
	2. If in the SRR is present _1 and _2 the sequencing is bidirectional
	3. Check if there are the same SRR repeated multiple times -> if yes in the point 2 they will be merged
	4. It will download a bash script to download the data
	5. Upload the script to the server
	```
	# Out of the server
	file='ena-file-download-read_run-PRJNA642125-fastq_ftp-20250317-1505.sh' # modify it according to your dataset
	scp /mnt/c/Users/carbi/Downloads/$file hbp694@med-bart.rd.tuni.fi:/data/omics/public_data/endocrine_disruptors/$GSE/scripts/
	#Put the password
	```
	6. Run the script
	```
	# In the server 
	tmux
	
	cd /data/omics/public_data/endocrine_disruptors/$GSE/raw_data/
	bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/$file
	
	#to enter in the conrol mode: ctrl+b 
	#to exit: d 
	#tmux ls
	#tmux a 
	```

## 2 - Convert SRR to GSM and Merging of the same GSM
To make easier the next mapping with the metadata, the SRR IDs are converted to the GSM (that will be found in a column of the metadata). 
```
# One time only:
#Download edirect or run the next 2 lines
cd
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo "export PATH=/home/hbp694/edirect:${PATH}" >> ${HOME}/.bashrc
export PATH=${HOME}/edirect:${PATH}
```

Run the script 
```
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/SRR_GSM.sh $GSE

# Quality chack before the merging
#conda activate seqtools
#cd /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/
#fastqc /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*
# Multiqc to combine multiple fastq reports in one overall report
#multiqc . 
```

If there are multiple equal GSM they will be merged in a unique one. 
```
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/Merge_GSM.sh $GSE
```

## 3 - Trimming, Alignment Unique mapping
Set up the environment with all the tools needed for the RNASeq pre-processing 
```
# One time only:
source /data/miniforge3/bin/activate
conda init 
#log out and log in again
```
**NOTE: now in front of the prompt should appeat the written "(base)"

Activate the environment
```
conda activate seqtools
```
**NOTE: before the prompt should appear the written "(seqtools)"

Run the quality check on the raw data with FASTQC.
```
cd /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/
fastqc /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/*
# Multiqc to combine multiple fastq reports in one overall report
multiqc . 
```
Check the fastqc and the multi fasqc to understand if the trimming is necessary and eventually which parameters to use. In the RNA-seq data, the Sequence Duplication Level is usually red because each transcript is present more than one time and also because after the PCR some replicates could be present. If the sequencing is performed with an Illumina machine it is possible that also the per base sequence content is red. These are not big deals at this stage if the other parameters are normal. 
The parameters to take into consideration are:
* Sequence Length (basic statistics) -> to understand the limit to put in the trimming part: usually the cut off is 60, but if the reads are shorter you can put 40. Not less than 40 otherwise the reads will be mapped in more than one spot in the genome.
* Per base sequence quality -> to understand how many bases have a quality score lower than 20 (usual threshold) 

```
#Change the parameters in the script according the information retrieved from the fastqc reports before running it
bash /data/omics/public_data/endocrine_disruptors/$GSE/scripts/trimming_alignment.sh $GSE
```

After the trimming with cutadapt a second quality control is performed. Check these FASTQCs even if the script will continue in any case. Then the alignment is performed against the Homo_sapiens.GRCh38_Ensembl_release_110 with hisat2. Finally the BAM files are processed to filter out all the reads with multiple mappings: we want to retain only the ones with one match against the genome. 

Look at the results of the alignment. If the percentage is low means that the alignment is not good. You have to ask yourself why.

## 4 - Counts Matrix
Change the conda environment
```
# Deactivate the old environment
conda deactivate
# Activate another environment to use R
conda activate renv
```
Run the script to have the counts Matrix
```
Rscript /data/omics/public_data/endocrine_disruptors/$GSE/scripts/count_matrix.R $GSE
```

** NOTE: 
1. Remember to delete the intermediate files and the raw data to free the space on the servers:
```
In all the directories:
rm *.bam
rm *.fastq.gz
rm *.fq.gz

#or

cd ..
find . -type f \( -name "*.bam" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) -exec rm {} \;
```
2. Use always 'tmux' so that no unexpected crash will happen
3. In bash remember to not put the space before and after the = when you declare the variables, otherwise it raises an error

##5 - Filtering and DEGs analysis
Script **filtering_analysis.Rmd**
Not automatized because you have to look at the single dataset. You have to look at the metadata, the controls and the plot to eventually remove the outliers 