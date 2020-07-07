# Nextflow RNA-seq pipeline
Nextflow pipeline to handle RNA-seq data (STAR 2-pass alignment, QC, featureCounts...)

https://github.com/maressyl/nextflow.RNA-seq



## Quick start


### Dependencies

* Singularity (tested with version `2.6.1-dist`)


### Building the Singularity container

`sudo singularity build RNA-seq.smig RNA-seq.def`


### Annotation files to download manually

- **Reference genome**, as a single FASTA file :
   - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
   
- **Reference transcriptome**, as a single GTF file :
   - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz

- **COSMIC coding mutations**, as a **bgzip-recompressed and indexed** VCF file :
   - https://cancer.sanger.ac.uk/cosmic/download

- **gnomAD polymorphisms**, as a bgzip compressed and indexed VCF file :
   - ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz


### Processing a run (example data provided)

```
#!/bin/bash

# Run parameters
FASTQ="data/test"

# Annotation files
FASTA="store/GRCh38.primary_assembly.genome.fa"
gnomAD="store/af-only-gnomad.hg38.vcf.gz"
COSMIC="store/CosmicCodingMuts_v90_GRCh38.vcf.gz"
GTF="store/gencode.v32.primary_assembly.annotation.gtf"

# Launch pipeline
nextflow -C "conf/base.conf" run "RNA-seq.nf" -with-singularity "RNA-seq.simg" --title "Test" --FASTQ "$FASTQ" \
  --readLength 76 --stranded 'R2' --RG_CN 'Integragen' --RG_PL 'ILLUMINA' --RG_PM 'HiSeq2000' \
  --genomeFASTA="$FASTA" --genomeGTF="$GTF" --gnomAD="$gnomAD" --COSMIC="$COSMIC" \
  --CPU_index 12 --CPU_align1 6 --CPU_align2 12 --CPU_mutect 12 \
  --out "./out" --single --varcall --window "chr7:148807000-148885000"
```
