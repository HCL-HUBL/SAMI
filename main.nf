#!/usr/bin/env nextflow

// Sample sheet (CSV, columns : sample, R1, R2)
params.input = ''
if(params.input == '') error "ERROR: --input must be provided"

// Series title
params.title = ''
if(params.title == '')                      error "ERROR: --title must be provided"
if(params.title ==~ /.*[^A-Za-z0-9_\.-].*/) error "ERROR: --title can only contain letters, digits, '.', '_' or '-'"

// Reference genome
params.species     = 'Human'
params.genome      = 'GRCh38'
params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
params.genomeFASTA = ''   // Example : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
params.genomeGTF   = ''   // Example : https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
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

// UMI-based read deduplication (optional)
params.umi = false
params.umi_protrude = 0

// SNV and indel calling (optional and experimental)
params.varcall = false
params.COSMIC = ''   // Example : https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v91/VCF/CosmicCodingMuts.vcf.gz + bgzip and .tbi index
params.gnomAD = ''   // Example : ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/af-only-gnomad.hg38.vcf.gz + bgzip and .tbi index
params.window = ''   // Genomic window into which restrict the variant calling (typically "chr7:148807000-148885000" to speed-up the test dataset)
if(params.varcall) {
	if(params.COSMIC == '') error "ERROR: --COSMIC must be provided with --varcall"
	if(params.gnomAD == '') error "ERROR: --gnomAD must be provided with --varcall"
	if(params.PL == '')     error "ERROR: --PL must be provided with --varcall"
}

// Aberrant splicing analysis
params.splicing = true
params.min_PSI = 0.1                    // Minimum "Percentage Spliced In" for an aberrant junction to be retained (between 0 and 1)
params.min_I = 30                       // Minimum reads supporting an aberrant junction to be retained
params.min_reads_unknown = 10           // "Unknown" junctions without this amount of reads or more in at least one sample will be ignored (significantly reduces computing time)
params.plot = true                      // Whether to plot genes with retained aberrant junctions or not
params.fusions = true                   // Whether to return gene fusions or ignore them
params.classes = "plausible,anchored"   // Classes of junctions to focus on during splicing analysis (comma-separated, among "unknown", "anchored", "plausible" and "annotated")
params.focus = "none"                   // IDs of junctions to focus on (chrom:start-end separated by commas), whatever their filtering status
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
	headerRegex = file("${projectDir}/modules/fastq/etc/FASTQ_headers.txt")
	fastq(
		FASTQ_pairs.groupTuple(by: 2),
		headerRegex,
		params.CN,
		params.PL,
		params.PM
	)

	// Build STAR index
	star_index(
		params.genomeFASTA,
		params.genomeGTF,
		params.genome
	)

	// STAR first pass
	star_pass1(
		fastq.out.FASTQ,
		star_index.out.genome,
		params.genomeGTF,
		params.umi_protrude
	)

	// Build a new genome from STAR pass 1
	dummy_R1 = file("${projectDir}/modules/STAR/reindex/etc/dummy_R1.fastq")
	dummy_R2 = file("${projectDir}/modules/STAR/reindex/etc/dummy_R2.fastq")
	star_reindex(
		star_pass1.out.junctions.collect(sort: true),
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
		params.genomeGTF,
		params.umi_protrude
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
		rrna_interval_genome.out.rRNA,
		stranded_Picard
	)
	rnaseqmetrics_target(
		bam_sort.out.BAM,
		"target",
		refflat_target.out.refFlat,
		rrna_interval_target.out.rRNA,
		stranded_Picard
	)

	// Count reads in transcripts using featureCounts
	featurecounts(
		bam_sort.out.BAM,
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
	secondary(bam_sort.out.BAM)

	// Plot soft-clipping lengths on read ends
	// TODO : general stats
	softclipping(bam_sort.out.BAM)

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
		fastqc_trimmed.out.ZIP.collect(sort: true),
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
		splicing_annotation(
			genomeGTF,
			params.species,
			params.genome,
			params.chromosomes
		)
		
		// Collect alignment gaps in each BAM
		splicing_harvest(
			bam_sort.out.BAM,
			indexfasta.out.indexedFASTA
		)
		
		// Transcript file channel (either used or empty file)
		if(params.transcripts != '') {
			transcripts = file(params.transcripts)
		} else {
			transcripts = file("${TMPDIR}/no-transcript.tsv")
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
			bam_sort.out.BAM.map{ it[2] }.collect(sort: true),
			bam_sort.out.BAM.map{ it[3] }.collect(sort: true),
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
		splicing_dir.add(params.focus.replaceAll(':','-'))
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
			params.classes,
			params.focus
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
			bqsr.out.BAM,
			params.window
		)
	}
}
