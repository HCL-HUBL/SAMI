#!/usr/bin/env bash

### YAML table
### Header
echo -e "custom_data:\n" \
    "    umi_duplication:\n" \
    "        plot_type: 'generalstats'\n" \
    "        pconfig:\n" \
    "            - UMI.duplication:\n" \
    "                namespace: 'UMI.duplication'\n" \
    "                description: 'Duplication based on the UMI (total read after consensus / total read before consensus)'\n" \
    "                format: '{:,.2f}'\n" \
    "        data:\n"  > duplication_umi.yaml

for read2 in ./*.consensus.fastq.gz
do
    sample=$(basename "${read2}" _R1_001.consensus.fastq.gz)
    read1=$(echo "${read2}" | sed 's/consensus\.//g' | sed 's/R1/L001_R1/g')
    ### Count all reads, even unmapped
    nread1=$(zcat "${read1}" | wc -l)
    nread2=$(zcat "${read2}" | wc -l)
    dup=$(echo "100-100*${nread2}/${nread1}" | bc -l)
    echo -e "            ${sample}:\n" \
        "                UMI.duplication: "${dup}""
done >> duplication_umi.yaml
