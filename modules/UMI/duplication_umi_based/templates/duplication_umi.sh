#!/usr/bin/env bash

### YAML table
### Header
echo -e "custom_data:\n" \
    "    umi_duplication:\n" \
    "        plot_type: 'generalstats'\n" \
    "        headers:\n" \
    "            - UMI.dup:\n" \
    "                namespace: 'UMI.dup'\n" \
    "                description: 'Proportion of reads considered as duplicates (1 - afterDedup / beforeDedup)'\n" \
    "                format: '{:,.2f}'\n" \
    "                min: 0\n" \
    "                max: 100\n" \
    "                suffix: '%'\n" \
    "        data:\n"  > duplication_umi.yaml

for read2 in ./*.DNA.MD.sort.bam
do
    sample=$(basename "${read2}" .DNA.MD.sort.bam)
    read1="./"$(basename "${read2}" .DNA.MD.sort.bam)".pass1.bam"
    ### Count without unmapped (0x0004), pair unmapped (0x0008) and secondary alignment (0x0100)
    nread1=$(samtools view -F 0x0004 -F 0x0008 -F 0x0100 -c "${read1}")
    nread2=$(samtools view -F 0x0004 -F 0x0008 -F 0x0100 -c "${read2}")
    dup=$(echo "100-100*${nread2}/${nread1}" | bc -l)

    echo -e "            ${sample}:\n" \
        "                UMI.dup: "${dup}""
done >> duplication_umi.yaml
