#!/usr/bin/env sh

INDIR=$1

### YAML table
### Header
echo -e "custom_data:\n" \
    "    duplication_umi:\n" \
    "        plot_type: 'generalstats'\n" \
    "        pconfig:\n" \
    "            - duplication.UMI:\n" \
    "                namespace: 'duplication.UMI'\n" \
    "                description: 'Duplication based on the UMI (nread after consensus / nread before consensus)'\n" \
    "                format: '{:,.2f}'\n" \
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
        "                duplication.UMI: "${dup}""
done >> duplication_umi.yaml
