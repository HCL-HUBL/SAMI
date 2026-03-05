#!/usr/bin/env nextflow

// Sample sheet (CSV, columns : sample, R1, R2)
params.input = ''
if(params.input == '') error "ERROR: --input must be provided"
params.fastq_check = true

// Series title
params.title = ''
if(params.title == '')                      error "ERROR: --title must be provided"
if(params.title ==~ /.*[^A-Za-z0-9_\.-].*/) error "ERROR: --title can only contain letters, digits, '.', '_' or '-'"

// Reference genome
params.species     = 'Human'
params.genome      = 'GRCh38'
params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
params.genomeFASTA = ''
params.genomeGTF   = ''
params.targetGTF   = ''
if(params.genomeFASTA == '') error "ERROR: --genomeFASTA must be provided"
if(params.genomeGTF == '')   error "ERROR: --genomeGTF must be provided"

// Read-group annotation (optional)
params.CN = ''
params.PL = 'ILLUMINA'
params.PM = ''

// Stranded library
params.stranded = 'no'
if(params.stranded == "R1") {
	stranded_Picard   = 'FIRST_READ_TRANSCRIPTION_STRAND'
	stranded_Rsubread = '1L'
} else if(params.stranded == "R2") {
	stranded_Picard   = 'SECOND_READ_TRANSCRIPTION_STRAND'
	stranded_Rsubread = '2L'
} else if(params.stranded == "no") {
	stranded_Picard   = 'NONE'
	stranded_Rsubread = '0L'
} else error "ERROR: --stranded must be 'R1', 'R2' or 'no'"

// Adapter trimming (optional)
params.trimR1 = ''
params.trimR2 = ''

// Adapter search with AdapteurRemoval (optional)
// If some triming values are provided and searchAdapteur
// is set to true, rise an error
params.identifyAdapter = false
if((params.trimR1 != '' || params.trimR2 != '') && params.identifyAdapter == true) error "ERROR: either --trimR1/--trimR2 or --identifyAdapter should be provided/set to true"

// UMI-based read deduplication (optional)
params.umi = false
params.umi_protrude = 0

// SNV and indel calling (optional and experimental)
params.varcall = false
params.gnomAD = ''
params.PoN = ''
params.vcf_format = 'full'
params.vcf_include = "FILTER='PASS'"
params.vcf_exclude = ''
if(params.varcall) {
	if(params.gnomAD == '') error "ERROR: --gnomAD must be provided with --varcall"
	if(params.PoN == '')    error "ERROR: --PoN must be provided with --varcall"
	if(params.PL == '')     error "ERROR: --PL must be provided with --varcall"
}

// Alignment
params.multimap = 5         // Maximum amount of mapping locations for a read to be considered aligned (-1 for all)
params.fixRange = 10        // Maximum distance to a known exon boundary to consider when trying to shift introns toward a single known splicing site
params.prefilter = 0.001    // Filter out during STAR pass 1 all junctions supported by less than this proportion of the sequencing depth at splicing site

// Aberrant splicing analysis
params.splicing = true
params.flags = 0                        // Any of these flags will exclude reads from junction counting (similar to samtools view -F ...)
params.min_PSI = 0.01                   // Minimum "Percentage Spliced In" for an aberrant junction to be retained (between 0 and 1)
params.min_I = 3                        // Minimum reads supporting an aberrant junction to be retained
params.min_reads_unknown = 10           // "Unknown" junctions without this amount of reads or more in at least one sample will be ignored (significantly reduces computing time)
params.plot = true                      // Whether to plot genes with retained aberrant junctions or not
params.fusions = true                   // Whether to return gene fusions or ignore them
params.classes = "plausible,anchored"   // Classes of junctions to focus on during splicing analysis (comma-separated, among "trivial", "nosplice", "unknown", "anchored", "plausible" and "annotated")
params.transcripts = ''                 // Preferred transcript table (2 tab-separated columns without header and quote : symbol and NCBI transcipt)
if(params.targetGTF == '') {
	// Symbols of genes to focus on during splicing analysis (comma-separated list, "all" to not filter or "target" to use symbols in targetGTF)
	params.symbols = "all"
} else {
	params.symbols = "target"
}

// Long-term storage
params.store = "./store"
params.output = "./output"

// Publish mode (how to deal with output files)
params.publish = "copy"

// Multi-QC annotation
params.MQC_title   = params.title
params.MQC_comment = ""



include { bcftools }                              from "./modules/bcftools"
include { cutadapt }                              from "./modules/adapter/cutadapt"
include { adapterremoval }                        from "./modules/adapter/adapterremoval"
include { retrieveadapter }                       from "./modules/adapter/adapterremoval"
include { fastq_check }                           from "./modules/fastq_check"
include { fastq_skip }                            from "./modules/fastq_skip"
include { featurecounts }                         from "./modules/featurecounts"
include { edgeR }                                 from "./modules/edgeR"
include { sample_sheet }                          from "./modules/sample_sheet"
include { star_index }                            from "./modules/STAR/index"
include { star_filtergaps }                       from "./modules/STAR/filtergaps"
include { star_fixgaps }                          from "./modules/STAR/fixgaps"
include { star_align as star_pass1 }              from "./modules/STAR/align"
include { star_align as star_pass2 }              from "./modules/STAR/align"
include { star_reindex }                          from "./modules/STAR/reindex"
include { indexfasta }                            from "./modules/Picard/indexfasta"
include { markduplicates }                        from "./modules/Picard/markduplicates"
include { bam_sort as sort_pass1 }                from "./modules/samtools/bam_sort"
include { bam_sort as sort_pass2 }                from "./modules/samtools/bam_sort"
include { depth as depth_pass1 }                  from "./modules/samtools/depth"
include { filterduplicates }                      from "./modules/samtools/filterduplicates"
include { umi_consensus }                         from "./modules/UMI/consensus"
include { duplication_umi_based }                 from "./modules/UMI/duplication_umi_based"
include { merge_filterbam }                       from "./modules/UMI/merge_filterbam"
include { umi_plot }                              from "./modules/UMI/plot"
include { umi_table }                             from "./modules/UMI/table"
include { bqsr }                                  from "./modules/GATK/bqsr"
include { indexvcf }                              from "./modules/GATK/indexvcf"
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
include { splicing_aggregate }                    from "./modules/splicing/aggregate"
include { splicing_annotation }                   from "./modules/splicing/annotation"
include { splicing_depth }                        from "./modules/splicing/depth"
include { splicing_filter }                       from "./modules/splicing/filter"
include { splicing_harvest }                      from "./modules/splicing/harvest"
include { splicing_nosplice }                     from "./modules/splicing/nosplice"



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

    if(params.identifyAdapter) {
		adapterremoval(FASTQ_pairs) // Identify the adapters for each pair
        retrieveadapter(adapterremoval.out.log.collect()) // Collect the log files and retrieve the adapters
        toTrimR1 = retrieveadapter.out.R1
        toTrimR2 = retrieveadapter.out.R2
    } else {
        // If no params.identifyAdapter, need to initialise toTrimR1/2
        toTrimR1 = params.trimR1
        toTrimR2 = params.trimR2
    }

	if(toTrimR1 != '' || toTrimR2 != '') {
        // Trim FASTQ
        // Use toTrimR1/2 to avoid initialising twice params.trimR1/2
		cutadapt(
			FASTQ_pairs,
			toTrimR1,
			toTrimR2
		)
		cutadapt_log = cutadapt.out.log.collect(sort: true)
		FASTQ_pairs = cutadapt.out.FASTQ
		
		// FastQC on trimmed FASTQ
		R1 = FASTQ_pairs.map{it[0]}.unique()
		R2 = FASTQ_pairs.filter{ it[4] == "paired" }.map{it[1]}.unique()
		fastqc_trimmed(
			R1.mix(R2)
		)
		fastqc_trimmed_ZIP = fastqc_trimmed.out.ZIP.collect(sort: true)
	} else {
		cutadapt_log = []
		fastqc_trimmed_ZIP = []
	}

	if(params.fastq_check) {
		// Check FASTQ headers and group by sample
		headerRegex = file("${projectDir}/modules/fastq_check/etc/FASTQ_headers.txt")
		fastq_check(
			FASTQ_pairs.groupTuple(by: 2),
			headerRegex,
			params.CN,
			params.PL,
			params.PM
		)
		FASTQ_pass1 = fastq_check.out.FASTQ
	} else {
		// Minimal FASTQ and group by sample
		fastq_skip(
			FASTQ_pairs.groupTuple(by: 2),
			params.CN,
			params.PL,
			params.PM
		)
		FASTQ_pass1 = fastq_skip.out.FASTQ
	}

	// Build STAR index
	star_index(
		params.genomeFASTA,
		params.genomeGTF,
		params.genome
	)
	
	// STAR first pass
	star_pass1(
		FASTQ_pass1,
		star_index.out.genome,
		params.genomeGTF,
		params.umi_protrude,
		params.multimap
	)
	
	if(params.prefilter > 0) {
		// Depth profile of STAR first pass
		sort_pass1(star_pass1.out.BAM_DNA)
		depth_pass1(sort_pass1.out.BAM)
		
		// Filter gaps from STAR first pass
		star_filtergaps(
			star_pass1.out.junctions.join(
				depth_pass1.out.TSV
			),
			params.prefilter
		)
		
		// Use filtered junctions
		junctions = star_filtergaps.out.junctions.collect(sort: true)
	} else {
		// Use raw junctions
		junctions = star_pass1.out.junctions.map{[ it[1] ]}.collect(sort: true)
	}
	
	// Prepare FASTA satellite files as requested by GATK
	indexfasta(params.genomeFASTA)
	
	// Prepare introns and exon track files
	splicing_annotation(
		params.genomeGTF,
		params.species,
		params.genome,
		params.chromosomes
	)
	
	// Collect and fix junctions from first pass
	star_fixgaps(
		splicing_annotation.out.exons,
		indexfasta.out.indexedFASTA,
		junctions,
		params.fixRange
	)
	
	// Build a new genome from STAR pass 1
	dummy_R1 = file("${projectDir}/modules/STAR/reindex/etc/dummy_R1.fastq")
	dummy_R2 = file("${projectDir}/modules/STAR/reindex/etc/dummy_R2.fastq")
	star_reindex(
		star_fixgaps.out.junctions,
		star_index.out.genome,
		params.genomeGTF,
		dummy_R1,
		dummy_R2,
		params.genome,
		params.title
	)

	if(params.umi) {
		// Create consensus reads from UMI-identified duplicates
		umi_consensus(
			star_pass1.out.BAM_DNA,
			params.CN,
			params.PL,
			params.PM
		)
		FASTQ_pass2 = umi_consensus.out.FASTQ
		
		// Convert duplication histogram for MultiQC
		umi_plot(
			umi_consensus.out.histogram
		)
		umi_plot_YAML = umi_plot.out.YAML.collect(sort: true)
		
		// Aggregate duplication table for MultiQC
		umi_table(
			umi_consensus.out.histogram.map{[ it[1] ]}.collect(sort: true)
		)
		umi_table_YAML = umi_table.out.YAML
	} else {
		// Use same reads as in pass 1
		FASTQ_pass2 = FASTQ_pass1
		
		umi_plot_YAML = []
		umi_table_YAML = []
	}

	// STAR second pass
	star_pass2(
		FASTQ_pass2,
		star_reindex.out.genome,
		params.genomeGTF,
		params.umi_protrude,
		params.multimap
	)
	
	// Estimate insert size distribution
	insertsize(star_pass2.out.isize)

	// Get the median insert size per sample
	insertsize_table(star_pass2.out.isize.filter { it[1] == "paired" }.map{it[2]}.collect(sort: true))

	if(params.umi) {
		// Merge and filter : consensus reads mapped + consensus reads unmapped + pass1 unmapped reads
		merge_filterbam(
			star_pass2.out.BAM_DNA.join(
				umi_consensus.out.BAM_unmapped.join(
					star_pass1.out.BAM_DNA.map{[ it[0], it[2] ]}
				)
			),
			indexfasta.out.indexedFASTA
		)
		BAM = merge_filterbam.out.BAM
	} else {
		// Use raw STAR BAM
		BAM = star_pass2.out.BAM_DNA
	}

	// Picard MarkDuplicates (mark only, filter later)
	// FIXME use as many CPUs as available, whatever the options
	// FIXME add a short @PG line (default adds to all reads and mess up with samtools other @PG)
	markduplicates(BAM)

	// Genomically sort and index
	sort_pass2(markduplicates.out.BAM)

	// Get duplication stats based on UMI
	if(params.umi) {
		duplication_umi_based(
			star_pass1.out.BAM_DNA.map{it[2]}.collect(sort: true),
			sort_pass2.out.BAM.map{it[2]}.collect(sort: true)
		)
		duplication_umi_based_YAML = duplication_umi_based.out.YAML
	} else {
		duplication_umi_based_YAML = []
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
		sort_pass2.out.BAM,
		"genome",
		refflat_genome.out.refFlat,
		rrna_interval_genome.out.rRNA,
		stranded_Picard
	)
	rnaseqmetrics_target(
		sort_pass2.out.BAM,
		"target",
		refflat_target.out.refFlat,
		rrna_interval_target.out.rRNA,
		stranded_Picard
	)

	// Count reads in transcripts using featureCounts
	featurecounts(
		sort_pass2.out.BAM,
		targetGTF,
		stranded_Rsubread
	)

	// Use edgeR to compute QC
	edgeR(
		featurecounts.out.annotation.first(),
		featurecounts.out.counts.collect(sort: true)
	)

	// Quantify secondary alignments with SAMtools
	// TODO : general stats
	secondary(sort_pass2.out.BAM)

	// Plot soft-clipping lengths on read ends
	// TODO : general stats
	softclipping(sort_pass2.out.BAM)

	// Collect QC files into a single report
	multiqc_conf = file("${projectDir}/modules/QC/multiqc/etc/multiqc.conf")
	multiqc(
		params.MQC_title,
		params.MQC_comment,
		multiqc_conf,
		edgeR.out.YAML_general,
		edgeR.out.YAML_section,
		star_pass1.out.log.collect(sort: true),
		star_pass2.out.log.collect(sort: true),
		fastqc_raw.out.ZIP.collect(sort: true),
		fastqc_trimmed_ZIP,
		markduplicates.out.txt.collect(sort: true),
		rnaseqmetrics_genome.out.RNA_Metrics.collect(sort: true),
		rnaseqmetrics_target.out.RNA_Metrics.collect(sort: true),
		insertsize.out.YAML.collect(sort: true),
		secondary.out.YAML.collect(sort: true),
		softclipping.out.YAML.collect(sort: true),
		umi_plot_YAML,
		umi_table_YAML,
		insertsize_table.out.YAML,
		cutadapt_log,
		duplication_umi_based_YAML,
		versions.out.YAML
	)

	if(params.splicing) {
		// Collect alignment gaps in each BAM
		splicing_harvest(
			sort_pass2.out.BAM,
			indexfasta.out.indexedFASTA,
			params.flags
		)
		
		// Transcript file channel (either used or empty file)
		if(params.transcripts != '') {
			transcripts = file(params.transcripts)
		} else {
			transcripts = file("/tmp/no-transcript.tsv")
			if(!transcripts.exists()) transcripts.text = ''
		}
		
		// Aggregate all splicing events
		splicing_aggregate(
			splicing_annotation.out.genes,
			splicing_annotation.out.exons,
			splicing_annotation.out.introns,
			splicing_harvest.out.TSV.collect(sort: true),
			star_pass2.out.chimeric.collect(sort: true),
			transcripts,
			params.chromosomes,
			params.min_reads_unknown,
			params.stranded
		)
		
		// Collect positions-of-interest sequencing depth in each BAM
		splicing_depth(
			sort_pass2.out.BAM.map{ it[2] }.collect(sort: true),
			sort_pass2.out.BAM.map{ it[3] }.collect(sort: true),
			splicing_aggregate.out.BED,
			10,
			30
		)
		
		// Add no-splice as an alternative to gaps
		splicing_nosplice(
			splicing_aggregate.out.RDS,
			splicing_depth.out.BED
		)
		
		// Output directory for splicing_filter
		splicing_dir = []
		splicing_dir.add("I-${params.min_I}")
		splicing_dir.add("PSI-${params.min_PSI}")
		splicing_dir.add("${params.symbols.take(50)}(${params.symbols.split(',').size()})")
		splicing_dir.add(params.classes)
		if(params.fusions) { splicing_dir.add("fusions")
		} else             { splicing_dir.add("no-fusions")
		}

		// Collect all splicing events
		splicing_filter(
			splicing_annotation.out.exons,
			splicing_nosplice.out.RDS,
			splicing_depth.out.BED,
			targetGTF,
			splicing_dir.join("_"),
			params.plot,
			params.fusions,
			params.min_I,
			params.min_PSI,
			params.symbols,
			params.classes
		)
	}

	// EXPERIMENTAL
	if(params.varcall) {
		if(params.umi) {
			// Deduplicate with UMIs
			varcall_BAM = sort_pass2.out.BAM
		} else {
			// Deduplicate with MarkDuplicates
			filterduplicates(sort_pass2.out.BAM)
			varcall_BAM = filterduplicates.out.BAM
		}

		// Picard SplitNCigarReads (split reads with intron gaps into separate reads)
		splitn(
			indexfasta.out.indexedFASTA,
			varcall_BAM
		)
		
		// Download and index VCF required by Mutect2
		indexvcf(
			params.gnomAD,
			params.PoN
		)
		
		// Call variants with GATK Mutect2
		mutect2(
			indexfasta.out.indexedFASTA,
			indexvcf.out.germline,
			indexvcf.out.PoN,
			splitn.out.BAM
		)
		
		// Convert VCF to TSV
		bcftools(
			mutect2.out.filtered_VCF,
			params.vcf_format,
			params.vcf_include,
			params.vcf_exclude
		)
	}
}
