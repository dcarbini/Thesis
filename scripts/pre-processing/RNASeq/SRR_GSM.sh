#!/bin/bash

##Convert SRR to GSM 

# Check if the user insert the argument 
if [ -z "$1" ]; then
    echo "Usage: $0 <GSE_ID>"
    exit 1
fi

# Assign the input to the variable
GSE="$1"

export EDIRECT_RETRY=10

# Raw_data directory
RAW_DIR="/data/omics/public_data/endocrine_disruptors/$GSE/raw_data"
# Output directory 
OUT_DIR="/data/omics/public_data/endocrine_disruptors/$GSE/raw_data_GSM"

# Determine file pattern (_1 or no _1)
if ls "$RAW_DIR"/*_1.fastq.gz >/dev/null 2>&1; then
    Files=$(ls "$RAW_DIR"/*_1.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $1}')
    paired=true
else
    Files=$(ls "$RAW_DIR"/*.fastq.gz | awk -F "/" '{print $NF}' | awk -F "[_.]" '{print $1}')
    paired=false
fi

# Initialize a counter
counter=1

for i in $Files; do
    # Extract the necessary parts from the experiment information
    experiment_info=$(esummary -db sra -id "$i" | grep -oP '<Title>\K[^<]+')

    # Extract the necessary parts from the experiment information
    experiment_info_short=$(echo "$experiment_info" | grep -oP 'GSM\d+')

    # Determine the number prefix based on the counter to avoid that equal GSM are overwritten
    prefix=$counter

    # Construct the new filenames
    if [ "$paired" = true ]; then
        new_filename_1="${prefix}_${experiment_info_short}_1.fastq.gz"
        new_filename_2="${prefix}_${experiment_info_short}_2.fastq.gz"

        cp "$RAW_DIR/${i}_1.fastq.gz" "$OUT_DIR/${new_filename_1}"
        cp "$RAW_DIR/${i}_2.fastq.gz" "$OUT_DIR/${new_filename_2}"
        
        echo "${i}_1 -> $new_filename_1"
        echo "${i}_2 -> $new_filename_2"
    else
        new_filename="${prefix}_${experiment_info_short}.fastq.gz"

        cp "$RAW_DIR/${i}.fastq.gz" "$OUT_DIR/${new_filename}"
        
        echo "${i} -> $new_filename"
    fi

    # Increment the counter for the next set of experiments
    counter=$((counter + 1))
done
