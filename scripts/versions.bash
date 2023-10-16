#!/bin/bash

echo "STAR: '$(STAR --version)'"
echo "Picard tools: '$(java -jar "$picard" MarkDuplicates --version 2>&1 | sed "s/Version://")'"
echo "fgbio: '$(java -jar "$fgbio" --version 2>&1 | sed -E "s/^.+Version: ([0-9\\.]+).+\$/\1/")'"
echo "EdgeR: '$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")'"
echo "samtools: '$(samtools --version | grep "Using htslib" | sed "s/Using //")'"
echo "Cutadapt: '$(cutadapt --version)'"
echo "FastQC: '$(fastqc -v | sed "s/FastQC v//")'"

