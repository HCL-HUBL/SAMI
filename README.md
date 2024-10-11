# Splicing Analysis with Molecular Indexes (SAMI)
A nextflow pipeline to handle RNA-seq data from FASTQ files to end results, with a special focus on splicing and gene-fusion events.

https://gitlab.inria.fr/HCL/pipelines/SAMI.git


![Overview](doc/SAMI_short.png)


## Table of content

[[_TOC_]]

## Quick start


### Dependencies

* Nextflow (tested with version `22.10.8-all`)
* Singularity (tested with version `CE 3.8.0`)


### Building the Singularity container

`sudo singularity build SAMI.sif SAMI.def`


### Annotation files to download (and gunzip) manually

#### GENCODE annotation

- Browse https://www.gencodegenes.org/human/ for latest versions.
- **Reference genome**, as a single FASTA file :
   - [GRCh38.primary_assembly.genome.fa.gz](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz)
- **Reference transcriptome**, as a single GTF file :
   - [gencode.v44.primary_assembly.annotation.gtf.gz](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz)

#### NCBI annotation

- Browse https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml for latest versions.
- **Reference genome**, as a single FASTA file :
   - [GCA_000001405.15_GRCh38_full_analysis_set.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz)
- **Reference transcriptome**, as a single GTF file :
   - [GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz)

#### Common files for variant-calling

- **COSMIC coding mutations**, as a **bgzip-recompressed and indexed** VCF file :
   - [CosmicCodingMuts.normal.vcf.gz](https://cancer.sanger.ac.uk/cosmic/archive-download) ("VCF files", registration required)

- **gnomAD polymorphisms**, as a bgzip compressed and indexed VCF file :
   - [af-only-gnomad.hg38.vcf.gz](https://www.bcgsc.ca/downloads/morinlab/reference/af-only-gnomad.hg38.vcf.gz)


### Processing a run (example data provided)

```bash
#!/bin/bash

# Annotation files (to be downloaded manually first)
genomeFASTA="store/GCA_000001405.15_GRCh38_full_analysis_set.fna"
genomeGTF="store/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"

# Launch pipeline
nextflow run SAMI.nf -with-singularity "SAMI.sif" \
   --FASTQ "./data" --out "./out" --title "SeraSeq" \
   --RG_CN 'HCL' --RG_PM 'MiSeq' --stranded "R2" \
   --genomeFASTA "$genomeFASTA" --genomeGTF "$genomeGTF"  \
   --CPU_index 12 --CPU_align1 6 --CPU_align2 12 \
   --CPU_cutadapt 2 --trimR1 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' --trimR2 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' \
   --umi --CPU_umi 6 --umi_length 6 \
   --splicing --CPU_splicing 12 --classes "plausible" --min_I 5 --min_PSI 0.1 --fusions true
```

## Argument description

### Main arguments

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--FASTQ | \<none but required\> | Path to a directory containing one sub-directory per sample, each containing a single `*_R1_001.fastq.gz` file and (optionnaly) a single `*_R2_001.fastq.gz` file. |
| \--title | \<none but required\> | Series’ name, to print on MultiQC report and use in file names. |
| \--genomeFASTA | \<none but required\> | Multi-FASTA file containing all chromosomes of the genomic reeference, for alignment (see example above). |
| \--genomeGTF | \<none but required\> | GTF file describing genes, transcripts and exons for the whole transcriptome (see example above). |
| \--targetGTF | \--genomeGTF | GTF file describing genes, transcripts and exons for the captured genes of interest (can be the whole transcriptome). |
| \--CPU\_index | \<none\> | How many CPUs to use to (re)build STAR indexes (recommended: 10). |
| \--CPU\_align1 | \<none\> | How many CPUs to use for STAR first pass (recommended: 5+). |
| \--CPU\_align2 | \<none\> | How many CPUs to use for STAR second pass (recommended: 5+). |
| \--species | "Human" | Name of the sample species, for file annotation. |
| \--genome | "GRCh38" | Name of the assembly of reference genome used, for file annotation. |
| \--chromosomes | "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" | Ordered list of chromosomes to consider during the analysis. |
| \--RG_CN | \<none\> | Sequencing center name (to populate the CN field in @RG of BAM files) |
| \--RG_PL | "ILLUMINA" | Sequencing technology used (to populate the PL field in @RG of BAM files, see SAM file specification for allowed values) |
| \--RG_PM | \<none\> | Sequencer model name (to populate the PM field in @RG of BAM files) |
| \--stranded | "no" | Whether a stranded RNA-seq library was used or not ("no", "R1" or "R2"), mainly used during QC. |
| \--store | "./store" | Path to long term storage for processed annotation files, to speed-up consecutive launchs of the pipeline. |
| \--out | "./out" | Path to output directory, where files of interest are published. |
| \--publish | "link" | Publishing mode for output files (see Nextflow documentation). |
| \--scratch | false | Temporary storage for processes (see Nextflow documentation). |
| \--MQC\_title | \--title | Title for the MultiQC report. |
| \--MQC\_comment | \<none\> | Free comment to add at the beginning of the MultiQC report. |
| \--single | false | Whether to handle single-end data (R1 only) or consider missing R2 file as an error. |

### Advanced arguments and debugging

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--clean\_BAM | true | Whether to empty intermediate BAM files in the 'work' directory to save storage space or not. |
| \--finalize | true | Whether to run MultiQC and EdgeR or not (for partial runs only). |
| \--fixedTime | \<none\> | Replace process time limits with this value (typically use '24h' to get rid of time limits). |
| \--maxRetries | 4 | Amount of times retriable processes with dynamic resources limits will restart when they encounter errors. |
| \--RNA\_BAM | false | Whether to publish the BAM file aligned to the transcriptome or not. |

### Adapter trimming (optional)

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--CPU\_cutadapt | \<none\> | How many CPUs to use for trimming (recommended: 2+). |
| \--trimR1 | \<none\> | Sequence to trim in 3’ of R1 (cutadapt -a). |
| \--trimR2 | \<none\> | Sequence to trim in 3’ of R2 (cutadapt -A). |

### UMI-based deduplication (optional)

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--umi | false | Whether to deduplicate reads based on pass 1 STAR alignment and UMI content (consensus read) or not. |
| \--CPU\_umi | \<none\> | How many CPUs to use for UMI processing (recommended: 6+). |
| \--umi\_length | \<none\> | Length of UMIs, only if they were located in 5’ of both R1 and R2 and extracted from the reads prior to launching SAMI (alignment parameters will be adjusted accordingly). Otherwise use 0. |

### Aberrant splicing analysis (optional)

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--splicing | true | Whether to look for aberrant splicing events or not. |
| \--CPU\_splicing | \<none\> | How many CPUs to use for splicing analysis (recommended: 5+). |
| \--min\_PSI | 0.1 | Minimum Percentage Spliced In (PSI) to retain an aberrant junction as a candidate (between 0 and 1). |
| \--min\_I | 30 | Minimum amount of (deduplicated) reads supporting an aberrant junction to retain it as a candidate. |
| \--min\_reads\_unknown | 10 | "Unknown" junctions without this amount of reads or more in at least one sample will be ignored (significantly reduces computing time). |
| \--plot | true | Whether to produce plots of genes and samples with at least one junction passing filters or not. |
| \--fusions | true | Whether to call gene-fusions or only splicing events inside genes. |
| \--classes | plausible | Classes of events to retain as candidates (see dedicated section below). |
| \--focus | "none" | Retain only these junctions as candidates (comma-separated list). |
| \--transcripts | \<none\> | Preferred transcripts to focus on (tab-separated file without header, with gene symbols in first column and transcript ID in second column). |
| \--symbols | "target" | Only junctions in these genes will be retained as candidates (use "all" to disable the filtering, or "target" to refer to genes defined in –targetGTF). |

### SNV and short-indel calling (experimental and optional)

| Argument | Default value | Description |
| :-- | :-- | :-- |
| \--varcall | false | Whether to perform SNV and short indel calling or not. |
| \--CPU\_mutect | \<none\> | How many CPUs to use for Mutect2 variant-calling (recommended: 4+). |
| \--COSMIC | \<none but required\> | VCF file of known pathogenic variants (bgzipped and TBI indexed) |
| \--gnomAD | \<none but required\> | VCF file of known polymorphisms (bgzipped and TBI indexed) |
| \--window | \<none\> | Genomic window in which to perform the variant calling (to speed-up tests mainly, leave empty to call in the entire genome). |

## UMI pre-processing

SAMI assumes Unique Molecular indexes (UMIs) are **not** part of the reads sequences, but were extracted (and concatenated) previously. Typically a read with a R1 UMI `CTTACT` and a R2 UMI `ATGTGA` would look like :

```
@M05299:437:000000000-K7V8N:1:1101:19520:2095:CTTACTATGTGA 1:N:0:GTGACAGT+CACTGTCA
CCCATATCTGCATTGTGACTAGCCTGTCATCCACCATCACCTCCTTGGTCAGAAAGTCAGCTCCTATTGTGGCTTTGTACTGATTGCTGAATTTCTTATTCACATACTGGTTCATGAGTGATGTCACATAGATCGGAAGAGCAC
+
@@>1D1FDG33AFGFGB3AFA1AAFCHG3GG11AB0FF1A0FEAFGF1FFF01111D2211FAAF1FGGGF?/FGFBGFFHH2GGHHHF11FGHFHHFHHFFFFGBFHB1FG11GB2BGF1GHHHHGHHB11GFC//?E00CFF
```

The easiest way to obtain such a result is to extract UMIs directly during demultiplexing, defining `U` cycles in `--use-bases-mask` (`bcl2fastq` software) or `OverrideCycles` (`bcl-convert` software).

To extract UMIs **a posteriori** in FASTQ files, software like [umi-tools](https://umi-tools.readthedocs.io/en/latest/) can be used. Be aware though that UMI-tools use `_` delimiter instead of `:` which need to be replaced manually to match SAMI's expectations.

## Events detected by SAMI

SAMI detects and classifies splicing events as follows :

### annotated

"annotated" junctions are introns described in at least one transcript of the provided annotation. They correspond to **physiological splicing**, useful to compute PSI values and assess gene coverage but not as candidate events.

### plausible

"plausible" junctions are splicing gaps joining the ends of two known exons, but this particular combination of exons is not found in any transcript of the provided annotation. They correspond typically to **exon skips**, and should be considered as high quality candidates.

### anchored

These gaps are "anchored" on one splicing site described in the provided annotation, either the left site (left-most genomic position, regardless of transcription strand) or the right site. The second splicing site involved in such an event is unknown to the annotation, and could correspond to a **neo-exon** or an **alternative splicing site**.

### unknown

Events are considered "unknown" if neither of the two involved splicing sites are described in the annotation. Such events are usually **alignment artefacts** but might be part of a **complex splicing event**. Notice they have to pass the extra `--min_reads_unknown` filter to be reported in the "All junctions" and "Candidates" files.