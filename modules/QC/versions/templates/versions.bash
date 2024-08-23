#!/bin/bash

out="SAMI_mqc_versions.yaml"
rm -f "$out"

# Nextflow variables
echo "Command: !{workflow.commandLine}"   >> "$out"
echo "Container: '!{workflow.container}'" >> "$out"
echo "Nextflow: '!{nextflow.version}'"    >> "$out"
echo "SAMI: '!{gitVersion}'"              >> "$out"

# Singularity-contained software
STAR --version                                                                 | xargs printf "STAR: '%s'\n" >> "$out"
java -jar "$picard" MarkDuplicates --version 2>&1 | sed "s/Version://"         | xargs printf "Picard tools: '%s'\n" >> "$out"
java -jar "$fgbio" --version 2>&1 | sed -E "s/^.+Version: ([0-9\\.]+).+\$/\1/" | xargs printf "fgbio: '%s'\n" >> "$out"
Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))"        | xargs printf "EdgeR: '%s'\n" >> "$out"
samtools --version | grep "Using htslib" | sed "s/Using //"                    | xargs printf "samtools: '%s'\n" >> "$out"
cutadapt --version                                                             | xargs printf "Cutadapt: '%s'\n" >> "$out"
fastqc -v | sed "s/FastQC v//"                                                 | xargs printf "FastQC: '%s'\n" >> "$out"
