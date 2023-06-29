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
    "                min: 0\n" \
    "                max: 100\n" \
    "                suffix: '%'\n" \
    "        data:\n"  > duplication_umi.yaml

for file in *_reads.txt
do
	sample=$(basename "${file%_reads.txt}")
	
	# Parse counts
	R1_before=$(grep "R1_before" $file | sed -E 's/^.+\t//')
	R1_after=$(grep "R1_after" $file | sed -E 's/^.+\t//')
	
	# Compute duplication
	dup=$(echo "100-100*${R1_after}/${R1_before}" | bc -l)
	
	# YAML content
	echo -e "            ${sample}:\n" \
		"                UMI.duplication: "${dup}""
done >> duplication_umi.yaml

