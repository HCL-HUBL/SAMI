#!/bin/bash

for sample in Sample_700611 Sample_700613 Sample_700650 Sample_700652
do
   echo $sample
   mkdir -p data/test/${sample}
   
   # 10 000 random reads
   for file in data/PRIMA/${sample}/*
   do
      gunzip --stdout $file | head -40000 | gzip --stdout > data/test/${sample}/$(basename $file)
   done
   
   # EZH2 locus
   samtools view -h "out.PRIMA/BAM/${sample}.DNA.sorted.bam" chr7:148807000-148885000 | samtools sort -n - | samtools fixmate - - | samtools fastq -N -f 0x2 -F 0x100 -F 0x800 -1 "data/test/${sample}/${sample}_EZH2_R1.fastq" -2 "data/test/${sample}/${sample}_EZH2_R2.fastq" -
   gzip data/test/${sample}/*_EZH2_*
done

# Make two samples single-end
for sample in Sample_700611 Sample_700652
do
   rm data/test/${sample}/*R2*
done

# Make a single-end sample single-file
gunzip data/test/Sample_700611/*
cat data/test/Sample_700611/*fastq | gzip --stdout > data/test/Sample_700611/Sample_700611_merged.R1.fq.gz
rm data/test/Sample_700611/*fastq

