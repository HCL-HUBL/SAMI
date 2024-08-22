#!/usr/bin/env nextflow

// Run characteristics (no default value)
params.input    = ''
params.stranded = ''
params.RG_CN    = ''
params.RG_PL    = 'ILLUMINA'
params.RG_PM    = ''
params.title    = ''

// Reference genome (files not provided)
params.species     = 'Human'
params.genome      = 'GRCh38'
params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
params.genomeFASTA = '' /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz */
params.genomeGTF   = '' /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz */
params.targetGTF   = '' /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz */
params.COSMIC      = '' /* https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v91/VCF/CosmicCodingMuts.vcf.gz + authentication / bgzip FIXME */
params.gnomAD      = '' /* ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz */

// CPU to use (no default value)
params.CPU_index    = 0
params.CPU_align1   = 0
params.CPU_align2   = 0
params.CPU_mutect   = 0
params.CPU_splicing = 0
params.CPU_cutadapt = 0
params.CPU_umi      = 0

// Whether to run splicing analysis or not
params.splicing = true

// Whether to run variant-calling processes or not
params.varcall = false

// FASTQ trimming sequences R1 and R2
params.trimR1 = ''
params.trimR2 = ''

// Need to generate UMI consensus read?
params.umi        = false
params.umi_length = -1

// Mandatory values (general)
if(params.input == '')                      error "ERROR: --input must be provided"
if(params.genomeFASTA == '')                error "ERROR: --genomeFASTA must be provided"
if(params.genomeGTF == '')                  error "ERROR: --genomeGTF must be provided"
if(params.CPU_index <= 0)                   error "ERROR: --CPU_index must be a positive integer (suggested: all available CPUs)"
if(params.CPU_align1 <= 0)                  error "ERROR: --CPU_align1 must be a positive integer (suggested: 6+)"
if(params.CPU_align2 <= 0)                  error "ERROR: --CPU_align2 must be a positive integer (suggested: 6+)"
if(params.title == '')                      error "ERROR: --title must be provided"
if(params.title ==~ /.*[^A-Za-z0-9_\.-].*/) error "ERROR: --title can only contain letters, digits, '.', '_' or '-'"

// Mandatory values (conditionnal)
if(params.umi) {
	if(params.CPU_umi <= 0)                 error "ERROR: --CPU_umi must be a positive integer (suggested: 6+) with --umi"
	if(params.umi_length < 0)               error "ERROR: --umi_length must be provided with --umi"
}
if(params.splicing) {
	if(params.CPU_splicing <= 0)            error "ERROR: --CPU_splicing must be a positive integer (suggested: 5+) with --splicing"
}
if(params.trimR1 != '' || params.trimR2 != '') {
	if(params.CPU_cutadapt <= 0)            error "ERROR: --CPU_cutadapt must be a positive integer (suggested: 2+) with --trimR1 or --trimR2"
}
if(params.varcall) {
	if(params.COSMIC == '')                 error "ERROR: --COSMIC must be provided with --varcall"
	if(params.gnomAD == '')                 error "ERROR: --gnomAD must be provided with --varcall"
	if(params.CPU_mutect <= 0)              error "ERROR: --CPU_mutect must be a positive integer (suggested: 4+) with --varcall"
	if(params.RG_PL == '')                  error "ERROR: --RG_PL must be provided with --varcall"
}

// Strandness
if(params.stranded == "R1") {
	params.stranded_Picard   = 'FIRST_READ_TRANSCRIPTION_STRAND'
	params.stranded_Rsubread = '1L'
} else if(params.stranded == "R2") {
	params.stranded_Picard   = 'SECOND_READ_TRANSCRIPTION_STRAND'
	params.stranded_Rsubread = '2L'
} else if(params.stranded == "no") {
	params.stranded_Picard   = 'NONE'
	params.stranded_Rsubread = '0L'
} else error "ERROR: --stranded must be 'R1', 'R2' or 'no'"

// Long-term storage
params.store = "${projectDir}/store"
params.out   = "${projectDir}/out"

// How to deal with output files (hard links by default, to safely remove the 'work' directory)
params.publish = "copy"

// Multi-QC annotation
params.MQC_title   = params.title
params.MQC_comment = ""

// To disable tailored per-process time limits, define a common time limit (typically '24h')
params.fixedTime = ''

// Maximum retry attempts for retriable processes with dynamic ressource limits
params.maxRetries = 4

// Genomic window into which restrict the variant calling (typically "chr7:148807000-148885000" to speed-up the test dataset)
params.window = ''

// Minimum "Percentage Spliced In" for an aberrant junction to be retained (between 0 and 1)
params.min_PSI = 0.1

// Minimum reads supporting an aberrant junction to be retained
params.min_I = 30

// "Unknown" junctions without this amount of reads or more in at least one sample will be ignored (significantly reduces computing time)
params.min_reads_unknown = 10

// Whether to plot genes with retained aberrant junctions or not
params.plot = true

// Whether to return gene fusions or ignore them
params.fusions = true

// Symbols of genes to focus on during splicing analysis (comma-separated list, "all" to not filter or "target" to use symbols in targetGTF)
if(params.targetGTF == '') {
	params.symbols = "all"
} else {
	params.symbols = "target"
}

// Classes of junctions to focus on during splicing analysis (comma-separated, among "unknown", "anchored", "plausible" and "annotated")
params.classes = "plausible"

// IDs of junctions to focus on (chrom:start-end separated by commas), whatever their filtering status
params.focus = "none"

// Preferred transcript table (2 tab-separated columns without header and quote : symbol and NCBI transcipt)
params.transcripts = ''

include { cutadapt }                              from "./modules/cutadapt"
include { fastq }                                 from "./modules/fastq"
include { featurecounts }                         from "./modules/featurecounts"
include { edgeR }                                 from "./modules/edgeR"
include { sample_sheet }                          from "./modules/sample_sheet"
include { star_index }                            from "./modules/STAR/index"
include { star_pass1 }                            from "./modules/STAR/pass1"
include { star_pass2 }                            from "./modules/STAR/pass2"
include { star_reindex }                          from "./modules/STAR/reindex"
include { indexfasta }                            from "./modules/Picard/indexfasta"
include { markduplicates }                        from "./modules/Picard/markduplicates"
include { bam_sort }                              from "./modules/samtools/bam_sort"
include { filterduplicates }                      from "./modules/samtools/filterduplicates"
include { umi_consensus }                         from "./modules/UMI/consensus"
include { duplication_umi_based }                 from "./modules/UMI/duplication_umi_based"
include { merge_filterbam }                       from "./modules/UMI/merge_filterbam"
include { umi_plot }                              from "./modules/UMI/plot"
include { umi_table }                             from "./modules/UMI/table"
include { bqsr }                                  from "./modules/GATK/bqsr"
include { mutect2 }                               from "./modules/GATK/mutect2"
include { splitn }                                from "./modules/GATK/splitn"
include { insertsize }                            from "./modules/QC/insertsize"
include { insertsize_table }                      from "./modules/QC/insertsize_table"
include { fastqc as fastqc_raw }                  from "./modules/QC/fastqc"
include { fastqc as fastqc_trimmed }              from "./modules/QC/fastqc"
include { multiqc }                               from "./modules/QC/multiqc"
include { refflat as refflat_genome }             from "./modules/QC/rnaseqmetrics"
include { refflat as refflat_target }             from "./modules/QC/rnaseqmetrics"
include { rnaseqmetrics as rnaseqmetrics_genome } from "./modules/QC/rnaseqmetrics"
include { rnaseqmetrics as rnaseqmetrics_target } from "./modules/QC/rnaseqmetrics"
include { rrna_interval as rrna_interval_genome } from "./modules/QC/rnaseqmetrics"
include { rrna_interval as rrna_interval_target } from "./modules/QC/rnaseqmetrics"
include { secondary }                             from "./modules/QC/secondary"
include { softclipping }                          from "./modules/QC/softclipping"
include { versions }                              from "./modules/QC/versions"
include { annotation }                            from "./modules/splicing/annotation"
include { splicing_collect }                      from "./modules/splicing/collect"
include { splicing_filter }                       from "./modules/splicing/filter"

workflow {
	// Collect software versions for MultiQC
	gitVersion = "git --git-dir=${projectDir}/.git describe --tags --long".execute().text.replaceAll("\\s","")
	versions(gitVersion)
	
	// FASTQ pair channel from sample sheet
	FASTQ_pairs = sample_sheet(params.input)
	
	// FastQC on raw FASTQ
	R1 = FASTQ_pairs.map{it[0]}.unique()
	R2 = FASTQ_pairs.filter{ it[4] == "paired" }.map{it[1]}.unique()
	fastqc_raw(
		R1.mix(R2)
	)
	
	if(params.trimR1 != '' || params.trimR2 != '') {
		// Trim FASTQ
		cutadapt(
			FASTQ_pairs,
			params.trimR1,
			params.trimR2
		)
		FASTQ_pairs = cutadapt.out.FASTQ
		
		// FastQC on trimmed FASTQ
		R1 = FASTQ_pairs.map{it[0]}.unique()
		R2 = FASTQ_pairs.filter{ it[4] == "paired" }.map{it[1]}.unique()
		fastqc_trimmed(
			R1.mix(R2)
		)
	}
	
	// Check FASTQ and group by sample
	headerRegex = Channel.value("$projectDir/in/FASTQ_headers.txt")
	fastq(
		FASTQ_pairs.groupTuple(by: 2),
		headerRegex
	)
	
	// Build STAR index
	star_index(
		params.genomeFASTA,
		params.genomeGTF
	)
	
	// STAR first pass
	star_pass1(
		fastq.out.FASTQ,
		star_index.out.genome,
		params.genomeGTF
	)
	
	// Build a new genome from STAR pass 1
	star_reindex(
		star_pass1.out.junctions.collect(sort: true),
		star_index.out.genome,
		params.genomeGTF
	)

	if(params.umi) {
		// Create consensus reads from UMI-identified duplicates
		umi_consensus(
			star_pass1.out.BAM_DNA
		)
		FASTQ_pass2 = umi_consensus.out.FASTQ
		
		// Convert duplication histogram for MultiQC
		umi_plot(
			umi_consensus.out.histogram
		)
		
		// Aggregate duplication table for MultiQC
		umi_table(
			umi_consensus.out.histogram.map{[ it[1] ]}.collect(sort: true)
		)
	} else {
		// Use same reads as in pass 1
		FASTQ_pass2 = fastq.out.FASTQ
	}
	
	// STAR second pass
	star_pass2(
		FASTQ_pass2,
		star_reindex.out.genome,
		params.genomeGTF
	)
	
	// Estimate insert size distribution
	insertsize(star_pass2.out.isize)

	// Get the median insert size per sample
	insertsize_table(star_pass2.out.isize.filter { it[1] == "paired" }.map{it[2]}.collect(sort: true))
	
	// Prepare FASTA satellite files as requested by GATK
	indexfasta(params.genomeFASTA)
	
	if(params.umi) {
		// Merge and filter : consensus reads mapped + consensus reads unmapped + pass1 unmapped reads
		merge_filterbam(
			star_pass2.out.BAM_DNA.join(
				umi_consensus.out.BAM_unmapped.join(
					star_pass1.out.BAM_DNA.map{[ it[0], it[1] ]}
				)
			),
			indexfasta.out.indexedFASTA
		)
		BAM = merge_filterbam.out.BAM
	} else {
		// Use STAR pass 2 BAM
		BAM = star_pass2.out.BAM
	}

	// Picard MarkDuplicates (mark only, filter later)
	// FIXME use as many CPUs as available, whatever the options
	// FIXME add a short @PG line (default adds to all reads and mess up with samtools other @PG)
	markduplicates(BAM)
	
	// Genomically sort and index
	bam_sort(markduplicates.out.BAM)
	
	// Get duplication stats based on UMI
	if(params.umi) {
		duplication_umi_based(
			star_pass1.out.BAM_DNA.map{it[1]}.collect(sort: true),
			bam_sort.out.BAM.map{it[2]}.collect(sort: true)
		)
	}
	
	// Prepare GTF files for preprocessing
	if(params.targetGTF == '') {
		genomeGTF = file(params.genomeGTF)
		targetGTF = file(params.genomeGTF)
	} else {
		genomeGTF = file(params.genomeGTF)
		targetGTF = file(params.targetGTF)
	}
	
	// Prepare refFlat file for Picard
	refflat_genome(genomeGTF)
	refflat_target(targetGTF)
	
	// Prepare rRNA interval list file for Picard
	rrna_interval_genome(
		genomeGTF,
		star_index.out.chrom
	)
	rrna_interval_target(
		targetGTF,
		star_index.out.chrom
	)
	
	// Picard's CollectRnaSeqMetrics
	rnaseqmetrics_genome(
		bam_sort.out.BAM,
		"genome",
		refflat_genome.out.refFlat,
		rrna_interval_genome.out.rRNA
	)
	rnaseqmetrics_target(
		bam_sort.out.BAM,
		"target",
		refflat_target.out.refFlat,
		rrna_interval_target.out.rRNA
	)
	
	// Count reads in transcripts using featureCounts
	featurecounts(
		bam_sort.out.BAM,
		targetGTF
	)

	// Use edgeR to compute QC
	edgeR(
		featurecounts.out.annotation.first(),
		featurecounts.out.counts.collect(sort: true)
	)
	
	// Quantify secondary alignments with SAMtools
	// TODO : general stats
	secondary(bam_sort.out.BAM)

	// Plot soft-clipping lengths on read ends
	// TODO : general stats
	softclipping(bam_sort.out.BAM)
	
	// Collect QC files into a single report
	multiqc(
		edgeR.out.YAML_general,
		edgeR.out.YAML_section,
		star_pass1.out.log.collect(sort: true),
		star_pass2.out.log.collect(sort: true),
		fastqc_raw.out.zip.collect(sort: true),
		fastqc_trimmed.out.zip.collect(sort: true),
		markduplicates.out.txt.collect(sort: true),
		rnaseqmetrics_genome.out.RNA_Metrics.collect(sort: true),
		rnaseqmetrics_target.out.RNA_Metrics.collect(sort: true),
		insertsize.out.YAML.collect(sort: true),
		secondary.out.YAML.collect(sort: true),
		softclipping.out.YAML.collect(sort: true),
		umi_plot.out.YAML.collect(sort: true),
		umi_table.out.YAML,
		insertsize_table.out.YAML,
		cutadapt.out.log.collect(sort: true),
		duplication_umi_based.out.YAML,
		versions.out.YAML
	)
	
	if(params.splicing) {
		// Prepare introns and exon track files
		annotation(genomeGTF)
		
		// Transcript file channel (either used or empty file)
		if(params.transcripts != '') {
			transcripts = Channel.value(params.transcripts)
		} else {
			transcripts = Channel.value("$projectDir/in/dummy.tsv")
		}
		
		// Collect all splicing events
		splicing_collect(
			annotation.out.genes,
			annotation.out.exons,
			annotation.out.introns,
			star_pass2.out.junctions.collect(sort: true),
			star_pass2.out.chimeric.collect(sort: true),
			transcripts
		)
		
		// Output directory for splicing_filter
		splicing_dir = []
		splicing_dir.add("I-${params.min_I}")
		splicing_dir.add("PSI-${params.min_PSI}")
		splicing_dir.add("${params.symbols.take(50)}(${params.symbols.split(',').size()})")
		splicing_dir.add(params.classes)
		splicing_dir.add(params.focus.replaceAll(':','-'))
		if(params.fusions) { splicing_dir.add("fusions")
		} else             { splicing_dir.add("no-fusions")
		}

		// Collect all splicing events
		splicing_filter(
			annotation.out.exons,
			splicing_collect.out.RDS,
			bam_sort.out.BAM.map{ it[2] }.collect(sort: true),
			bam_sort.out.BAM.map{ it[3] }.collect(sort: true),
			targetGTF,
			splicing_dir.join("_")
		)
	}
	
	// EXPERIMENTAL
	if(params.varcall) {
		// Filter out duplicated read, based on a previous MarkDuplicates run
		filterduplicates(bam_sort.out.BAM)
	
		// Picard SplitNCigarReads (split reads with intron gaps into separate reads)
		splitn(
			indexfasta.out.indexedFASTA,
			filterduplicates.out.BAM
		)
		
		// SNV references
		gnomAD = Channel.value( [ params.gnomAD , params.gnomAD + ".tbi" ] )
		COSMIC = Channel.value( [ params.COSMIC , params.COSMIC + ".tbi" ] )
		
		// Compute and apply GATK Base Quality Score Recalibration model
		bqsr(
			indexfasta.out.indexedFASTA,
			gnomAD,
			COSMIC,
			splitn.out.BAM
		)

		// Call variants with GATK Mutect2
		// FIXME : --panel-of-normals pon.vcf.gz
		mutect2(
			indexfasta.out.indexedFASTA,
			gnomAD,
			bqsr.out.BAM
		)
	}
}
