#!/usr/bin/env nextflow

// Run characteristics (no default value)
params.input    = ''
params.stranded = ''
params.RG_CN    = ''
params.RG_PL    = ''
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
params.splicing = false

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

// Temporary storage ('work' directory by default, process memory directives don't account for 'ram-disk' usage !)
params.scratch = "false"

// How to deal with output files (hard links by default, to safely remove the 'work' directory)
params.publish = "copy"

// Last git commit (for versioning)
gitVersion = "git --git-dir=${projectDir}/.git describe --tags --long".execute().text.replaceAll("\\s","")

// Multi-QC annotation
params.MQC_title   = params.title
params.MQC_comment = ""

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

// Include the processes
include { versions }               from "./modules/QC/versions"
include { fastq }                  from "./modules/fastq"
include { cutadapt }               from "./modules/cutadapt"
include { fastqc_raw }             from "./modules/QC/fastqc_raw"
include { fastqc_trimmed }         from "./modules/QC/fastqc_trimmed"
include { star_index }             from "./modules/STAR/index"
include { star_pass1 }             from "./modules/STAR/pass1"
include { umi_stat_and_consensus } from "./modules/UMI/stat_and_consensus"
include { umi_plot }               from "./modules/UMI/plot"
include { umi_table }              from "./modules/UMI/table"
include { star_reindex }           from "./modules/STAR/reindex"
include { star_pass2 }             from "./modules/STAR/pass2"
include { insertsize_table }       from "./modules/QC/insertsize_table"
include { indexfasta }             from "./modules/Picard/indexfasta"
include { merge_filterbam }        from "./modules/UMI/merge_filterbam"
include { markduplicates }         from "./modules/Picard/markduplicates"
include { bam_sort }               from "./modules/samtools/bam_sort"
include { duplication_umi_based }  from "./modules/UMI/duplication_umi_based"
include { filterduplicates }       from "./modules/samtools/filterduplicates"
include { splitn }                 from "./modules/GATK/splitn"
include { bqsr }                   from "./modules/GATK/bqsr"
include { mutect2 }                from "./modules/GATK/mutect2"
include { refflat }                from "./modules/QC/rnaseqmetrics"
include { rrna_interval }          from "./modules/QC/rnaseqmetrics"
include { rnaseqmetrics }          from "./modules/QC/rnaseqmetrics"
include { featurecounts }          from "./modules/featurecounts"
include { edgeR }                  from "./modules/edgeR"
include { insertsize }             from "./modules/QC/insertsize"
include { secondary }              from "./modules/QC/secondary"
include { softclipping }           from "./modules/QC/softclipping"
include { multiqc }                from "./modules/QC/multiqc"
include { annotation }             from "./modules/splicing/annotation"
include { splicing_collect }       from "./modules/splicing/collect"
include { splicing_filter }        from "./modules/splicing/filter"

workflow {
	// Collect R1 and R2 per sample
	FASTQ_map = [:];
	sampleSheet = file("${params.input}")
	lines = sampleSheet.splitCsv(header: true)
	for(line in lines) {
		sampleName = line["sample"]
		if(FASTQ_map.containsKey(sampleName)) {
			FASTQ_map[sampleName]["R1"] << line["R1"]
			FASTQ_map[sampleName]["R2"] << line["R2"]
		} else {
			FASTQ_map[sampleName] = [ "R1": [ line["R1"] ], "R2": [ line["R2"] ] ]
		}
	}
	
	// Pairing type
	FASTQ_list = [];
	FASTQ_map.each { sampleName, sample ->
		// Check lists
		if(sample["R1"].size() != sample["R2"].size()) {
			error "ERROR: R1 and R2 file counts differ for $sampleName"
		}
		if(sample["R1"].size() == 0) {
			error "ERROR: No R1 file for $sampleName"
		}
		
		// Check elements
		anyPE = false
		for(int i = 0; i < sample["R1"].size(); i++) {
			if(sample["R1"][i] == "") {
				error "ERROR: Empty R1 file path for $sampleName"
			}
			if(sample["R2"][i] == "") {
				if(anyPE) {
					error "ERROR: Mixed single and paired end files for $sampleName"
				}
			} else {
				anyPE = true
			}
		}
		
		// Paired-end or single-end
		if(anyPE) { type = "paired"
		} else    { type = "single"
		}
		
		// Reformat as a list
		FASTQ_list << [ "R1": sample["R1"], "R2": sample["R2"], "sample": sampleName, "type": type ]
	}
	
	// Reformat as a channel
	FASTQ = Channel.fromList(FASTQ_list)
	
	// No insertSize output is OK (only single-end data)
	insertSize_bypass = Channel.fromPath("${projectDir}/in/dummy.tsv")

	// Annotation file channels
	if(params.targetGTF == '') {
		targetGTF = Channel.value(params.genomeGTF)
	} else {
		targetGTF = Channel.value(params.targetGTF)
	}
	genomeGTF   = Channel.value(params.genomeGTF)
	genomeFASTA = Channel.value(params.genomeFASTA)
	headerRegex = Channel.value("$projectDir/in/FASTQ_headers.txt")
	if(params.varcall) {
		gnomAD_Mutect2 = Channel.value( [ params.gnomAD , params.gnomAD + ".tbi" ] )
		COSMIC         = Channel.value( [ params.COSMIC , params.COSMIC + ".tbi" ] )
	} else {
		gnomAD_Mutect2 = Channel.of()
		COSMIC         = Channel.of()
	}

	// Transcript file channel (either used or empty file)
	if(params.splicing && params.transcripts != '') {
		transcripts = Channel.value(params.transcripts)
	} else {
		transcripts = Channel.value("$projectDir/in/dummy.tsv")
	}

	// Collect software versions for MultiQC
	versions(gitVersion)

	// Build RG line from 1st read of each FASTQ file pair bundle
	fastq(FASTQ,
		  headerRegex)

	// Run cutadapt
	// or bypass
	if(params.trimR1 != '' || params.trimR2 != '') {
		cutadapt(fastq.out.FASTQ_CUTADAPT)
	} else {
		cutadapt.out.outR1       = Channel.fromPath("${R1}")
		cutadapt.out.outR2       = Channel.fromPath("${R2}")
		cutadapt.out.outQC       = Channel.fromPath("${projectDir}/in/dummy.tsv")
		cutadapt.out.FASTQ_STAR1 = fastq.out.FASTQ_CUTADAPT
	}

	// Run FastQC on individual FASTQ files (raw FASTQ)
	fastqc_raw(fastq.out.R1_raw.concat(fastq.out.R2_raw))

	// Run FastQC on individual FASTQ files
	// or bypass FastQC_trimmed
	if(params.trimR1 != '' || params.trimR2 != '') {
		fastqc_trimmed(cutadapt.out.R1_trimmed.concat(
			cutadapt.out.R2_trimmed))
	} else {
		fastqc_trimmed.out.outQC = Channel.fromPath("${projectDir}/in/dummy.tsv")
	}

	// Build STAR index
	star_index(genomeFASTA,
			   genomeGTF)

	// STAR first pass
	// TODO shared-memory
	star_pass1(cutadapt.out.FASTQ_STAR1,
			   star_index.out.rawGenome,
			   genomeGTF)

	// If UMI present, take care of them
	// otherwise bypass
	if(params.umi) {
		// Change read name, the "_" into a ":" before the UMI in read name;
		// Create an unmapped BAM and mapped it with STAR
		// see: https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md
		umi_stat_and_consensus(star_pass1.out.BAM_pass1,
							   cutadapt.out.FASTQ_STAR1)

		// Get the UMI duplication stat in the FASTQC
		umi_plot(umi_stat_and_consensus.out.UMI_stat)

		// Get the UMI table
		umi_table(umi_stat_and_consensus.out.UMI_histo.collect())
	} else {
		umi_stat_and_consensus.out.FASTQ_STAR2 = star_pass1.out.FASTQ_STAR1_copy
		umi_plot.out.QC_umi                    = Channel.fromPath("${projectDir}/in/dummy.tsv")
		umi_table.out.QC_umi_table             = Channel.fromPath("${projectDir}/in/dummy.tsv")
	}

	// Build a new genome from STAR pass 1
	star_reindex(star_pass1.out.SJ_pass1.collect(),
				 genomeGTF,
				 star_index.out.rawGenome)

	// STAR second pass
	// TODO shared-memory
	star_pass2(umi_stat_and_consensus.out.FASTQ_STAR2,
			   star_reindex.out.reindexedGenome,
			   genomeGTF)

	// Get the median insert size per sample
	insertsize_table(star_pass2.out.isize_table.collect())

	// Prepare FASTA satellite files as requested by GATK
	indexfasta(genomeFASTA)

	// Merge mapped and unmapped BAM and filter
	// or skip it if no UMI
	if(params.umi) {
		merge_filterbam(star_pass2.out.genomic_temp_BAM.join(
			umi_stat_and_consensus.out.BAM_unmapped.join(
				star_pass1.out.BAM_pass1)),
						indexfasta.out.indexedFASTA)
	} else {
		"""
		mv "${BAM_mapped}" "${sample}.DNA.bam"
		"""
		merge_filterbam.out.genomic_BAM = "${sample}.DNA.bam"
	}

	// Picard MarkDuplicates (mark only, filter later)
	// FIXME use as many CPUs as available, whatever the options
	// FIXME add a short @PG line (default adds to all reads and mess up with samtools other @PG)
	markduplicates(merge_filterbam.out.genomic_BAM)

	// Genomically sort and index
	bam_sort(markduplicates.out.BAM_marked)

	// Get duplication stats based on UMI
	// or bypass
	if(params.umi) {
		duplication_umi_based(star_pass1.out.BAM_dup1.collect(),
							  bam_sort.out.onlyBAM_sorted.collect())
	} else {
		duplication_umi_based.out.outYAML = Channel.fromPath("$projectDir/in/dummy.tsv")
	}

	// Filter out duplicated read, based on a previous MarkDuplicates run
	filterduplicates(bam_sort.out.BAM_sorted)

	// Picard SplitNCigarReads (split reads with intron gaps into separate reads)
	splitn(indexfasta.out.indexedFASTA,
		   filterduplicates.out.BAM_filtered)

	// Compute and apply GATK Base Quality Score Recalibration model
	bqsr(indexfasta.out.indexedFASTA,
		 gnomAD_Mutect2,
		 COSMIC,
		 splitn.out.BAM_splitN)

	// Call variants with GATK Mutect2 (FIXME : --panel-of-normals pon.vcf.gz)
	mutect2(indexfasta.out.indexedFASTA,
			gnomAD_Mutect2,
			bqsr.out.BAM_BQSR)

	// Prepare refFlat file for Picard
	refflat(targetGTF.mix(genomeGTF).unique())

	// Prepare rRNA interval list file for Picard
	rrna_interval(targetGTF.mix(genomeGTF)
				  .unique(),
				  star_index.out.rawGenome_chrom)

	// Picard's CollectRnaSeqMetrics
	rnaseqmetrics(bam_sort.out.BAM_sorted
				  .combine(rrna_interval.out.rRNAs)
				  .combine(refflat.out.refFlats),
				  genomeGTF,
				  targetGTF)

	// Count reads in transcripts using featureCounts
	featurecounts(bam_sort.out.BAM_sorted,
				  targetGTF)

	// Use edgeR to compute QC
	edgeR(featurecounts.out.featureCounts_annotation.first(),
		  featurecounts.out.featureCounts_counts.collect())

	// Estimate insert size distribution
	insertsize(star_pass2.out.isize_sample)

	// Quantify secondary alignments with SAMtools
	// TODO : general stats
	secondary(bam_sort.out.BAM_sorted)

	// Plot soft-clipping lengths on read ends
	// TODO : general stats
	softclipping(bam_sort.out.BAM_sorted)

	// Collect QC files into a single report
	multiqc(edgeR.out.QC_edgeR_general,
			edgeR.out.QC_edgeR_section,
			star_pass1.out.QC_STAR_pass1.collect(),
			star_pass2.out.QC_STAR_pass2.collect(),
			fastqc_raw.out.QC_FASTQC_raw.collect(),
			fastqc_trimmed.out.QC_FASTQC_trimmed.collect(),
			markduplicates.out.QC_markDuplicates.collect(),
			rnaseqmetrics.out.QC_rnaSeqMetrics.collect(),
			insertSize_bypass.mix(insertsize.out.QC_insert).collect(),
			secondary.out.QC_secondary.collect(),
			softclipping.out.QC_softClipping.collect(),
			umi_plot.out.QC_umi.collect(),
			umi_table.out.QC_umi_table,
			insertsize_table.out.median_isize_table,
			cutadapt.out.QC_cutadapt.collect(),
			duplication_umi_based.out.dup_umi,
			versions.out.versions)

	// Prepare introns and exon track files
	annotation(genomeGTF)

	// Collect all splicing events
	splicing_collect(annotation.out.genes,
					 annotation.out.exons,
					 annotation.out.introns,
					 star_pass2.out.junctions_STAR.collect(),
					 star_pass2.out.chimeric_STAR.collect(),
					 transcripts)

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
	splicing_filter(annotation.out.exons,
					splicing_collect.out.splicing_events,
					bam_sort.out.onlyBAM_sorted.collect(),
					bam_sort.out.BAI_splicing.collect(),
					targetGTF,
					splicing_dir.join("_"))
}
