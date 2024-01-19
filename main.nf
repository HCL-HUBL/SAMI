#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { version } from "./modules/version"
include { fastq } from "./modules/fastq"
include { cutadapt } from "./modules/cutadapt"
include { fastqc_raw } from "./modules/fastqc_raw"
include { fastqc_trimmed } from "./modules/fastqc_trimmed"
include { star_index } from "./modules/star_index"
include { star_pass1 } from "./modules/star_pass1"
include { umi_stat_and_consensus } from "./modules/umi_stat_and_consensus"
include { umi_plot } from "./modules/umi_plot"
include { umi_table } from "./modules/umi_table"
include { star_reindex } from "./modules/star_reindex"
include { star_pass2 } from "./modules/star_pass2"
include { insertsize_table } from "./modules/insertsize_table"
include { indexfasta } from "./modules/indexfasta"
include { merge_filterbam } from "./modules/merge_filterbam"
include { markduplicates } from "./modules/markduplicates"
include { bam_sort } from "./modules/bam_sort"
include { duplication_umi_based } from "./modules/duplication_umi_based"
include { filterduplicates } from "./modules/filterduplicates"
include { splitN } from "./modules/splitn"
include { bqsr } from "./modules/bqsr"
include { mutect2 } from "./modules/mutect2"
include { refflat } from "./modules/refflat"
include { rrna_interval } from "./modules/rrna_interval"
include { rnaseqmetrics } from "./modules/rnaseqmetrics"
include { featurecounts } from "./modules/featurecounts"
include { edgeR } from "./modules/edgeR"
include { insertsize } from "./modules/insertsize"
include { secondary } from "./modules/secondary"
include { softclipping } from "./modules/softclipping"
include { multiqc } from "./modules/multiqc"
include { clean_dna_bam } from "./modules/clean_dna_bam"
include { clean_rna_bam } from "./modules/clean_rna_bam"
include { annotation } from "./modules/annotation"
include { splicing_collect } from "./modules/splicing_collect"
include { splicing_filter } from "./modules/splicing_filter"

workflow {

	// Run characteristics (no default value)
	params.FASTQ = ''
	params.stranded = ''
	params.RG_CN = ''
	params.RG_PL = ''
	params.RG_PM = ''
	params.title = ''

	// Reference genome (files not provided)
	params.species = 'Human'
	params.genome = 'GRCh38'
	params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
	params.genomeFASTA = ''   /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz */
	params.genomeGTF = ''     /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz */
	params.targetGTF = ''     /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz */
	params.COSMIC = ''        /* https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v91/VCF/CosmicCodingMuts.vcf.gz + authentication / bgzip FIXME */
	params.gnomAD = ''        /* ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz */

	// CPU to use (no default value)
	params.CPU_index = 0
	params.CPU_align1 = 0
	params.CPU_align2 = 0
	params.CPU_mutect = 0
	params.CPU_splicing = 0
	params.CPU_cutadapt = 0
	params.CPU_umi = 0

	// Whether to run splicing analysis or not
	params.splicing = false

	// Whether to run variant-calling processes or not
	params.varcall = false

	// FASTQ trimming sequences R1 and R2
	params.trimR1 = ''
	params.trimR2 = ''

	// Need to generate UMI consensus read?
	params.umi = false
	params.umi_length = -1

	// Mandatory values (general)
	if(params.FASTQ == '')                      error "ERROR: --FASTQ must be provided"
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
	}

	// Strandness
	if(params.stranded == "R1") {
		params.stranded_Picard = 'FIRST_READ_TRANSCRIPTION_STRAND'
		params.stranded_Rsubread = '1L'
	} else if(params.stranded == "R2") {
		params.stranded_Picard = 'SECOND_READ_TRANSCRIPTION_STRAND'
		params.stranded_Rsubread = '2L'
	} else if(params.stranded == "no") {
		params.stranded_Picard = 'NONE'
		params.stranded_Rsubread = '0L'
	} else error "ERROR: --stranded must be 'R1', 'R2' or 'no'"



	// Long-term storage
	params.store = "${baseDir}/store"
	params.out = "${baseDir}/out"

	// Temporary storage ('work' directory by default, process memory directives don't account for 'ram-disk' usage !)
	params.scratch = 'false'

	// How to deal with output files (hard links by default, to safely remove the 'work' directory)
	params.publish = 'link'

	// Last git commit (for versioning)
	gitVersion = "git --git-dir=${baseDir}/.git describe --tags --long".execute().text.replaceAll("\\s","")

	// Multi-QC annotation
	params.MQC_title = params.title
	params.MQC_comment = ""

	// Whether to publish BAM files aligning to the transcriptome or not
	params.RNA_BAM = false

	// Whether to remove unnecessary BAM files (unpublished RNA.bam and intermediary DNA.bam) from work or not (experimental)
	params.clean_BAM = true

	// To enable final processes assuming all samples were included (MultiQC and edgeR)
	params.finalize = true

	// To disable tailored per-process time limits, define a common time limit (typically '24h')
	params.fixedTime = ''

	// Maximum retry attempts for retriable processes with dynamic ressource limits
	params.maxRetries = 4

	// Whether to handle single-end data (R1 only) or consider missing R2 file as an error
	params.single = false

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

	// Collect FASTQ files from sample-specific folders
	FASTQ_list = []
	fastqDirectory = path("${params.FASTQ}")
	fastqDirectory.eachDir { sampleDirectory ->
		sample = sampleDirectory.name
		R1 = []
		R2 = []

		// For each R1 file
		anySE = false
		anyPE = false
		sampleDirectory.eachFileMatch(~/.*_R1_001\.fastq.gz/) { R1_file ->
			// FIXME add arguments for more flexibility (R1/R3 and pattern)
			// Corresponding R3 file (if any, assumes it is R2)
			R3_name = R1_file.name.replaceFirst(/(.*)_R1_001\.fastq.gz/, '$1_R3_001.fastq.gz')
			R3_file = file("${params.FASTQ}/${sample}/${R3_name}")
			if(R3_file.exists()) {
				// Use R3 as R2
				R2_file = R3_file;
				anyPE = true
			} else {
				// Corresponding R2 file
				R2_name = R1_file.name.replaceFirst(/(.*)_R1_001\.fastq.gz/, '$1_R2_001.fastq.gz')
				R2_file = path("${params.FASTQ}/${sample}/${R2_name}")
				if(R2_file.exists()) {
					// Use R2 as R2
					anyPE = true
				} else if(params.single) {
					// Neither R2 nor R3 : consider as single-end
					R2_file = ""
					anySE = true
				} else {
					// Neither R2 nor R3 : consider as an error
					error "ERROR: missing R2 file '${R2_name}' for sample '${sample}' (no R3 neither)"
				}
			}

			// Collect files
			R1.add(R1_file)
			R2.add(R2_file)
		}

		// Single or paired ends
		if(anyPE && anySE) error "ERROR: mixed single-end and paired-end samples ('${sample}') are not handled"
		if(anyPE) type = "paired"
		if(anySE) type = "single"

		// "Empty" directory
		if(R1.size() == 0) error "ERROR: no R1 file detected for sample '${sample}'"

		// Send to the channel
		FASTQ_list << [ "R1": R1, "R2": R2, "sample": sample, "type": type ]
	}
	FASTQ = Channel.fromList(FASTQ_list)

	// No insertSize output is OK (only single-end data)
	insertSize_bypass = Channel.of('dummy')

	// Annotation file channels
	if(params.targetGTF == '') {
		targetGTF = Channel.value(path(params.genomeGTF))
	} else {
		targetGTF = Channel.value(path(params.targetGTF))
	}
	genomeGTF = Channel.value(path(params.genomeGTF))
	genomeFASTA = Channel.value(path(params.genomeFASTA))
	headerRegex = Channel.value(path("$baseDir/in/FASTQ_headers.txt"))
	if(params.varcall) {
		gnomAD_BQSR    = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
		gnomAD_Mutect2 = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
		gnomAD_Mutect2 = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
		COSMIC         = Channel.value( [ file(params.COSMIC) , file(params.COSMIC + ".tbi") ] )
		//	rawCOSMIC      = Channel.value(file(params.COSMIC))
	} else {
		gnomAD_BQSR    = Channel.of()
		gnomAD_Mutect2 = Channel.of()
		COSMIC         = Channel.of()
		//	rawCOSMIC      = Channel.of()
	}

	// Transcript file channel (either used or empty file)
	if(params.splicing && params.transcripts != '') {
		transcripts = Channel.value(path(params.transcripts))
	} else {
		transcripts = Channel.value(path("$baseDir/in/dummy.tsv"))
	}

	// Collect software versions for MultiQC
	version()

	// Build RG line from 1st read of each FASTQ file pair bundle
	fastq(FASTQ,
		  headerRegex)

	// Run cutadapt
	cutadapt(fastq.FASTQ_CUTADAPT)

	// Run FastQC on individual FASTQ files (raw FASTQ)
	fastqc_raw(fastq.R1_raw.concat(fastq.R2_raw))

	// Run FastQC on individual FASTQ files
	fastqc_trimmed(fastqc_cutadapt
				   .R1_trimmed
				   .concat(fastqc_trimmed
						   .R2_trimmed))

	// Build STAR index
	star_index(genomeFASTA,
			   genomeGTF)

	// STAR first pass
	// TODO shared-memory
	star_pass1(cutadapt.FASTQ_STAR1,
			   star_index.rawGenome_pass1,
			   genomeGTF)

	// Change read name, the "_" into a ":" before the UMI in read name;
	// Create an unmapped BAM and mapped it with STAR
	// see: https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md
	umi_stat_and_consensus(star_pass1.BAM_pass1,
						   star_pass1.FASTQ_STAR1_copy)

	// Get the UMI duplication stat in the FASTQC
	umi_plot(umi_stat_and_consensus.UMI_stat)

	// Get the UMI table
	umi_table(umi_stat_and_consensus.UMI_table.collect())

	// Build a new genome from STAR pass 1
	star_reindex(star_pass1.SJ_pass1.collect(),
				 genomeGTF,
				 star_index.rawGenome_reindex)

	// STAR second pass
	// TODO shared-memory
	star_pass2(umi_stat_and_consensus.FASTQ_STAR2,
			   star_reindex.reindexedGenome,
			   genomeGTF)

	// Get the median insert size per sample
	insertsize_table(star_pass2.isize_table.collect())

	// Prepare FASTA satellite files as requested by GATK
	indexfasta(genomeFASTA)

	// Merge mapped and unmapped BAM and filter
	// or skip it if no UMI
	merge_filterbam(star_pass2.genomic_temp_BAM.join(umi_stat_and_consensus.BAM_unmapped.join(star_pass1.BAM_forUnmappedRead)),
					index_fasta.indexedFASTA_MergeBamAlignment)

	// Picard MarkDuplicates (mark only, filter later)
	// FIXME use as many CPUs as available, whatever the options
	// FIXME add a short @PG line (default adds to all reads and mess up with samtools other @PG)
	markduplicates(merge_filterbam.genomic_BAM)

	// Genomically sort and index
	bam_sort(markduplicates.BAM_marked)

	// Get duplication stats based on UMI
	duplication_umi_based(star_pass1.BAM_dup1.collect(),
						  bam_sort.BAM_dup2.collect())

	// Filter out duplicated read, based on a previous MarkDuplicates run
	filterduplicates(bam_sort.BAM_sorted)

	// Picard SplitNCigarReads (split reads with intron gaps into separate reads)
	splitn(indexfasta.indexedFASTA_splitN,
		   filterduplicates.BAM_filtered)

	// Compute and apply GATK Base Quality Score Recalibration model
	bqsr(indexfasta.indexedFASTA_BQSR,
		 gnomAD_BQSR,
		 COSMIC,
		 splitn.BAM_splitN)

	// Call variants with GATK Mutect2 (FIXME : --panel-of-normals pon.vcf.gz)
	mutect2(indexfasta.indexedFASTA_Mutect2,
			gnomAD_Mutect2,
			bqsr.BAM_BQSR)

	// Prepare refFlat file for Picard
	refflat(targetGTF.mix(genomeGTF).unique())

	// Prepare rRNA interval list file for Picard
	rrna_interval(targetGTF.mix(genomeGTF)
				  .unique(),
				  star_index.rawGenome_chrom)

	// Picard's CollectRnaSeqMetrics
	rnaseqmetrics(bam_sort.BAM_rnaSeqMetrics
				  .combine(rrna_interval.rRNAs)
				  .combine(refflat.refFlats))

	// Count reads in transcripts using featureCounts
	featurecounts(bam_sort.BAM_featureCounts,
				  targetGTF)

	// Use edgeR to compute QC
	edgeR(featurecounts.featureCounts_annotation.first(),
		  featurecounts.featureCounts_counts.collect())

	// Estimate insert size distribution
	insertsize(star_pass2.isize_sample)

	// Quantify secondary alignments with SAMtools
	// TODO : general stats
	secondary(bam_sort.BAM_secondary)

	// Plot soft-clipping lengths on read ends
	// TODO : general stats
	softclipping(bam_sort.BAM_softClipping)

	// Collect QC files into a single report
	multiqc(edgeR.QC_edgeR_general,
			edgeR.QC_edgeR_section,
			star_pass1.QC_STAR_pass1.collect(),
			star_pass2.QC_STAR_pass2.collect(),
			fastqc_raw.QC_FASTQC_raw.collect(),
			fastqc_trimmed.QC_FASTQC_trimmed.collect(),
			markduplicates.QC_markDuplicates.collect(),
			rnaseqmetrics.QC_rnaSeqMetrics.collect(),
			insertSize_bypass.mix(insertsize.QC_insert).collect(),
			secondary.QC_secondary.collect(),
			softclipping.QC_softClipping.collect(),
			umi_plot.QC_umi.collect(),
			umi_table.QC_umi_table,
			insertsize_table.median_isize_table,
			cutadapt.QC_cutadapt.collect(),
			duplication_umi_based.dup_umi,
			versions.versions)

	// Remove unnecessary BAM (unstorable process)
	clean_dna_bam(bam_sort.BAM_sort_clean
				  .mix(markduplicates.markDuplicates_clean,
					   splitn.splitN_clean,
					   bqsr.BQSR_clean))

	// Remove unnecessary BAM (unstorable process)
	clean_rna_bam(star_pass2.transcriptomic_BAM)

	// Prepare introns and exon track files
	annotation(genomeGTF)

	// Collect all splicing events
	splicing_collect(genes,
					 annotation.exons_collect,
					 annotation.introns,
					 star_pass2.junctions_STAR.collect(),
					 star_pass2.chimeric_STAR.collect(),
					 transcriptsl)

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
	splicing_filter(annotation.exons_filter,
					splicing_collect.splicing_events,
					bam_sort.BAM_splicing.collect(),
					bam_sort.BAI_splicing.collect(),
					targetGTF,
					splicing_dir.join("_"))
}
