#!/bin/bash

## Merge the same GSM 

# Check if the user insert the argument 
if [ -z "$1" ]; then
    echo "Usage: $0 <GSE_ID>"
    exit 1
fi

# Assign the input to the variable
GSE="$1"

# Get unique GSMS values
GSM_Values=$(ls /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*GSM*.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $2}' | sort | uniq)

for GSM_Value in $GSM_Values; do
    # Check if the files have _1 and _2
    if ls /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*"$GSM_Value"_1.fastq.gz >/dev/null 2>&1; then
        cat /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*"$GSM_Value"_1.fastq.gz > /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/"$GSM_Value"_1.fastq.gz
        cat /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*"$GSM_Value"_2.fastq.gz > /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/"$GSM_Value"_2.fastq.gz
    else
        cat /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM/*"$GSM_Value".fastq.gz > /data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM_merged/"$GSM_Value".fastq.gz
    fi
done
