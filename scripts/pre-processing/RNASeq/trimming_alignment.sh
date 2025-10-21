#!/bin/bash


# Check if the user insert the argument 
if [ -z "$1" ]; then
    echo "Usage: $0 <GSE_ID>"
    exit 1
fi

# Assign the input to the variable
GSE="$1"

# Determine file pattern (_1 or no _1) and get the list of FASTQ files
if ls /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/*_1.fastq.gz >/dev/null 2>&1; then
    FASTQ_FILES=$(ls /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/*_1.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $1}')
    paired=true
else
    FASTQ_FILES=$(ls /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/*.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $1}')
    paired=false
fi

# Run cutadapt and save logs
for i in $FASTQ_FILES; do
	if [ "$paired" = true ]; then
		cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 -m 60 -j 40 \
		-o /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}_1.fastq.gz \
		-p /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}_2.fastq.gz \
		/data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/${i}_1.fastq.gz \
		/data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/${i}_2.fastq.gz | tee /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/${i}_cutadapt.log
	else
		cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -q 20 -m 60 -j 40 \
		-o /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}.fastq.gz \
		/data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/${i}.fastq.gz | \
		tee /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/${i}_cutadapt.log

	fi
done
#Parameters to adapt: adapters (-a, -A) if the sequencing machine is not an Illumina one; length of the reads (-m); used cores (-j)

echo "Cutadapt Finished"


# Run FastQC on trimmed data
cd /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/
fastqc /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_* 
multiqc .

echo "FASTQC Finished"

# Get list of trimmed FASTQ files
FASTQ_FILES=$(ls /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_*.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $2}')


# Run hisat2 and save logs
for i in $FASTQ_FILES; do
	if [ "$paired" = true ]; then
		hisat2 -q -p 30 -x /data/omics/references/hsapiens/primary_assembly/hisat2_indexes/Homo_sapiens.GRCh38_Ensembl_release_110 \
	   -1 /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}_1.fastq.gz \
	   -2 /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}_2.fastq.gz \
		2> >(tee /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/${i}_hisat2.log >&2) | \
	   samtools view -Sbh > /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/TRIMMED_${i}.bam
	else
	   hisat2 -q -p 30 -x /data/omics/references/hsapiens/primary_assembly/hisat2_indexes/Homo_sapiens.GRCh38_Ensembl_release_110 \
	   -U /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/trimmed_fastq/TRIMMED_${i}.fastq.gz \
		2> >(tee /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/${i}_hisat2.log >&2) | \
	   samtools view -Sbh > /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/TRIMMED_${i}.bam
	fi
done

echo "Alignment Finished"

# Process BAM files
BAMFILES=$(ls /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/TRIMMED*.bam | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $2}')

for i in $BAMFILES; do
    samtools view -H /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/TRIMMED_${i}.bam > /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_header.sam
    samtools view /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/bam/TRIMMED_${i}.bam | grep -w "NH:i:1" | cat /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_header.sam - | \
    samtools view -@ 10 -Sb - > /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_unique.bam
    samtools sort -@ 10 /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_unique.bam > /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_unique_sorted.bam
    rm /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_header.sam
    rm /data/omics/public_data/endocrine_disruptors/$GSE/processed_data/unique_sorted/${i}_unique.bam
done

echo "BAM files processed"