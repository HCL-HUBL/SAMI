# Nextflow RNA-seq pipeline
Nextflow pipeline to handle RNA-seq data (STAR 2-pass alignment, QC, featureCounts...)

https://gitlab.inria.fr/NGS/pipelines/rna-seq.git


## Quick start


### Dependencies

* Nextflow (tested with version `21.12.1-edge-all`)
* Singularity (tested with version `CE 3.8.0`)


### Building the Singularity container

`sudo singularity build RNA-seq.sif RNA-seq.def`


### Annotation files to download (and gunzip) manually

#### GENCODE annotation

- Browse https://www.gencodegenes.org/human/ for latest versions.
- **Reference genome**, as a single FASTA file :
   - [GRCh38.primary_assembly.genome.fa.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz)
- **Reference transcriptome**, as a single GTF file :
   - [gencode.v41.primary_assembly.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz)

#### NCBI annotation

- Browse https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml for latest versions.
- **Reference genome**, as a single FASTA file :
   - [GCA_000001405.15_GRCh38_full_analysis_set.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz)
- **Reference transcriptome**, as a single GTF file :
   - [GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz)

#### Common files for variant-calling

- **COSMIC coding mutations**, as a **bgzip-recompressed and indexed** VCF file :
   - https://cancer.sanger.ac.uk/cosmic/download

- **gnomAD polymorphisms**, as a bgzip compressed and indexed VCF file :
   - ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz


### Processing a run (example data provided)

```bash
#!/bin/bash

# Run parameters
FASTQ="data/test"

# Annotation files
FASTA="store/GRCh38.primary_assembly.genome.fa"
gnomAD="store/af-only-gnomad.hg38.vcf.gz"
COSMIC="store/CosmicCodingMuts_v90_GRCh38.vcf.gz"
GTF="store/gencode.v32.primary_assembly.annotation.gtf"

# Launch pipeline
nextflow -C "conf/base.conf" run "RNA-seq.nf" -with-singularity "RNA-seq.sif" --title "Test" --FASTQ "$FASTQ" \
  --readLength 76 --stranded 'R2' --RG_CN 'Integragen' --RG_PL 'ILLUMINA' --RG_PM 'HiSeq2000' \
  --genome="GRCh38" --genomeFASTA="$FASTA" --genomeGTF="$GTF" --gnomAD="$gnomAD" --COSMIC="$COSMIC" \
  --CPU_index 12 --CPU_align1 6 --CPU_align2 12 --CPU_mutect 12 \
  --out "./out" --single --varcall --window "chr7:148807000-148885000" \
  --splicing --CPU_splicing 12 --refGene="store/refGene.txt.gz" --min_I=10 --min_PSI=0.05
```
