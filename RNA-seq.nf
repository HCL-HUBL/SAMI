#!/usr/bin/env nextflow

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
lastCommit = "git --git-dir=${baseDir}/.git log --format='%h' -n 1".execute().text.replaceAll("\\s","")

// Multi-QC annotation
params.MQC_title = params.title
params.MQC_comment = "Processed with maressyl/nextflow.RNA-seq [ ${lastCommit} ]"

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
fastqDirectory = file("${params.FASTQ}")
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
			R2_file = file("${params.FASTQ}/${sample}/${R2_name}")
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
FASTQ = Channel.from(FASTQ_list)

// No insertSize output is OK (only single-end data)
insertSize_bypass = Channel.from('dummy')

// Annotation file channels
if(params.targetGTF == '') {
	targetGTF = Channel.value(file(params.genomeGTF))
} else {
	targetGTF = Channel.value(file(params.targetGTF))
}
genomeGTF = Channel.value(file(params.genomeGTF))
genomeFASTA = Channel.value(file(params.genomeFASTA))
headerRegex = Channel.value(file("$baseDir/in/FASTQ_headers.txt"))
if(params.varcall) {
	gnomAD_BQSR    = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
	gnomAD_Mutect2 = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
	gnomAD_Mutect2 = Channel.value( [ file(params.gnomAD) , file(params.gnomAD + ".tbi") ] )
	COSMIC         = Channel.value( [ file(params.COSMIC) , file(params.COSMIC + ".tbi") ] )
//	rawCOSMIC      = Channel.value(file(params.COSMIC))
} else {
	gnomAD_BQSR    = Channel.from()
	gnomAD_Mutect2 = Channel.from()
	COSMIC         = Channel.from()
//	rawCOSMIC      = Channel.from()
}

// Transcript file channel (either used or empty file)
if(params.splicing && params.transcripts != '') {
	transcripts = Channel.value(file(params.transcripts))
} else {
	transcripts = Channel.value(file("$baseDir/in/dummy.tsv"))
}

// Build RG line from 1st read of each FASTQ file pair bundle
process FASTQ {
	
	cpus 1
	label 'nonRetriable'
	
	// Never scratch to avoid full copy of output in ram-disk
	scratch false
	stageInMode 'symlink'
	executor 'local'

	input:
	set file(R1), file(R2), val(sample), val(type) from FASTQ
	file regex from headerRegex
	
	output:
	set file(R1), file(R2), val(sample), val(type), stdout into FASTQ_CUTADAPT
	file R1 into R1_raw
	file R2 into R2_raw

	"""
	#!/usr/bin/env Rscript --vanilla
	
	# Get FASTQ sets from Nextflow (FIXME not space-proof)
	R1 <- strsplit("$R1", split=" ", fixed=TRUE)[[1]]
	R2 <- strsplit("$R2", split=" ", fixed=TRUE)[[1]]
	
	# Parse regex list
	def <- scan("${regex}", what="", sep="\n", quiet=TRUE)
	regex <- def[ c(TRUE,FALSE) ]
	names <- strsplit(def[ c(FALSE,TRUE) ], split=" ")
	
	# For each R1/R2 pair
	RG <- NULL
	for(i in 1:length(R1)) {
		
		# Single-end
		pairedEnd <- "${type}" == "paired"
		
		# Get first read headers (whether the file is compressed or not)
		H1 <- scan(R1[i], what="", sep="\n", n=1L, quiet=TRUE)
		H2 <- scan(R2[i], what="", sep="\n", n=1L, quiet=TRUE)
		if(length(H1) == 0L)              stop("No header in R1")
		if(length(H2) == 0L && pairedEnd) stop("No header in R2")
		
		# For each defined regex
		metadata_R1 <- character(0)
		metadata_R2 <- character(0)
		for(j in 1:length(regex)) {
			# Matching regex (R1)
			if(length(metadata_R1) == 0L && grepl(regex[j], H1)) {
				# Extract elements
				metadata_R1 <- regmatches(H1, regexec(regex[j], H1))[[1]][-1]
				names(metadata_R1) <- names[[j]]
			}
			
			# Matching regex (R2)
			if(pairedEnd && length(metadata_R2) == 0L && grepl(regex[j], H2)) {
				# Extract elements
				metadata_R2 <- regmatches(H2, regexec(regex[j], H2))[[1]][-1]
				names(metadata_R2) <- names[[j]]
			}
		}
		
		# No match
		if(length(metadata_R1) == 0L)              stop("Unable to parse the header of R1")
		if(pairedEnd && length(metadata_R2) == 0L) stop("Unable to parse the header of R2")
		
		# Run identity
		ID1 <- metadata_R1[ c("instrument", "run", "flowcell", "lane", "index") ]
		if(pairedEnd) {
			ID2 <- metadata_R2[ c("instrument", "run", "flowcell", "lane", "index") ]
			if(!identical(ID1, ID2)) stop("Wrong R1/R2 FASTQ file matching")
		}
		
		# Mate identity (if any is provided, both are required)
		R1_has_read_info <- !is.na(metadata_R1["read"]) && metadata_R1["read"] != ""
		R2_has_read_info <- !is.na(metadata_R2["read"]) && metadata_R2["read"] != ""
		if(R1_has_read_info || R2_has_read_info) {
			if(is.na(metadata_R1["read"]) || metadata_R1["read"] != "1")                   stop("R1 FASTQ file does not contain R1 reads")
			if(pairedEnd && (is.na(metadata_R2["read"]) || !metadata_R2["read"] %in% 2:3)) stop("R2 FASTQ file does not contain R2 or R3 reads")
		}
		
		# BC is optional
		if(is.na(metadata_R1["index"]) || metadata_R1["index"] == "") { BC <- ""
		} else                                                        { BC <- sprintf("\tBC:%s", metadata_R1["index"])
		}
		
		# RG definition for BAM header (current pair)
		x <- sprintf("ID:%s_%i%s", "${sample}", i, BC)
		if("${params.RG_CN}" != "") x <- c(x, "CN:${params.RG_CN}")
		if("${params.RG_PL}" != "") x <- c(x, "PL:${params.RG_PL}")
		if("${params.RG_PM}" != "") x <- c(x, "PM:${params.RG_PM}")
		x <- c(x, sprintf("PU:%s", paste(ID1, collapse=":")))
		x <- c(x, "SM:${sample}")
		
		# Merge with all read pairs
		RG <- c(RG, paste(x, collapse="\t"))
	}
	
	# Print final RG to stdout
	cat(paste(RG, collapse=" , "))
	"""
}

if(params.trimR1 != '' || params.trimR2 != '') {
	// Run cutadapt
	process cutadapt {

		cpus { params.CPU_cutadapt }
		label 'monocore'
		label 'retriable'

		storeDir { "${params.out}/cutadapt" }

		input:
		set file(R1), file(R2), val(sample), val(type), val(RG) from FASTQ_CUTADAPT

		output:
		set file("${R1.getSimpleName()}_cutadapt.fastq.gz"), file("${R2.getSimpleName()}_cutadapt.fastq.gz"), val(sample), val(type), val(RG) into FASTQ_STAR1
		file "${R1.getSimpleName()}_cutadapt.fastq.gz" into R1_trimmed
		file "${R2.getSimpleName()}_cutadapt.fastq.gz" into R2_trimmed
		file "${sample}_cutadapt.log" into QC_cutadapt

		"""
		if [[ ${params.trimR1} != "" ]] && [[ ${params.trimR2} != "" ]]
		then
			cutadapt -j ${params.CPU_cutadapt} \
				-a "${params.trimR1}" \
				-A "${params.trimR2}" \
				--minimum-length 20 \
				-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
				-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
				"${R1}" "${R2}" > "${sample}_cutadapt.log"
		elif [[ ${params.trimR1} != "" ]]
		then
			cutadapt -j ${params.CPU_cutadapt} \
				-a "${params.trimR1}" \
				--minimum-length 20 \
				-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
				-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
				"${R1}" "${R2}" > "${sample}_cutadapt.log"
		else
			cutadapt -j ${params.CPU_cutadapt} \
				-A "${params.trimR2}" \
				--minimum-length 20 \
				-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
				-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
				"${R1}" "${R2}" > "${sample}_cutadapt.log"
		fi
		"""
	}
} else {
	// Bypass cutadapt
	QC_cutadapt = Channel.value(file("$baseDir/in/dummy.tsv"))
	FASTQ_CUTADAPT.set{ FASTQ_STAR1 }
	R1_trimmed = Channel.from()
	R2_trimmed = Channel.from()
}

// Run FastQC on individual FASTQ files (raw FASTQ)
process FastQC_raw {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	
	storeDir { "${params.out}/QC/FastQC/raw" }
	
	when:
	!(FASTQ.name =~ /__input\.[0-9]+$/)
	
	input:
	file(FASTQ) from R1_raw.concat(R2_raw)
	
	output:
	file "${FASTQ.getSimpleName()}_fastqc.zip" into QC_FASTQC_raw
	
	"""
	fastqc "$FASTQ" -o "."
	"""
}

// Run FastQC on individual FASTQ files
if(params.trimR1 != '' || params.trimR2 != '') {
	process FastQC_trimmed {

		cpus 1
		label 'monocore'
		label 'retriable'

		storeDir { "${params.out}/QC/FastQC/trimmed" }

		when:
		!(FASTQ.name =~ /__input\.[0-9]+$/)

		input:
		file(FASTQ) from R1_trimmed.concat(R2_trimmed)

		output:
		file "${FASTQ.getSimpleName()}_fastqc.zip" into QC_FASTQC_trimmed

		"""
		fastqc "$FASTQ" -o "."
		"""
	}
} else {
	// Bypass FastQC_trimmed
	QC_FASTQC_trimmed = Channel.value(file("$baseDir/in/dummy.tsv"))
}

// Build STAR index
process STAR_index {
	
	cpus { params.CPU_index }
	label 'multicore'
	label 'nonRetriable'
	storeDir { params.store }
	
	input:
	file genomeFASTA from genomeFASTA
	file genomeGTF from genomeGTF
	
	output:
	file("${params.genome}_raw") into (rawGenome_pass1, rawGenome_reindex)
	file("${params.genome}_raw/chrNameLength.txt") into rawGenome_chrom
	
	"""
	mkdir -p "./${params.genome}_raw"
	STAR \
		--runThreadN ${params.CPU_index} \
		--runMode genomeGenerate \
		--genomeDir "./${params.genome}_raw" \
		--genomeFastaFiles "$genomeFASTA" \
		--sjdbGTFfile "$genomeGTF"
	"""
}

// STAR first pass
// TODO shared-memory
process STAR_pass1 {
	
	cpus { params.CPU_align1 }
	label 'multicore'
	label 'retriable'
	storeDir { "${params.out}/STAR_pass1" }
	
	input:
	set file(R1), file(R2), val(sample), val(type), val(RG) from FASTQ_STAR1
	file rawGenome from rawGenome_pass1
	file genomeGTF from genomeGTF

	output:
	file("${sample}.SJ.out.tab") into SJ_pass1
	set file("${sample}.pass1.bam"), val(sample), val(type), val(RG) into BAM_pass1
	set file(R1), file(R2), val(sample), val(type), val(RG) into FASTQ_STAR1_copy
	file "${sample}.pass1.bam" into BAM_dup1
	set val(sample), file("${sample}.pass1.bam") into BAM_forUnmappedRead
	file "${sample}_Log.final.out" into QC_STAR_pass1

	"""
	mkdir -p "./$sample"
	
	# FASTQ files
	if [ "$type" = "paired" ];   then readFilesIn="\\"${R1.join(",")}\\" \\"${R2.join(",")}\\""
	elif [ "$type" = "single" ]; then readFilesIn="\\"${R1.join(",")}\\""
	else                         echo "Unknow type '$type'"; exit 1
	fi
	
	STAR \
		--runThreadN ${params.CPU_align1} \
		--twopassMode None \
		--alignEndsProtrude ${params.umi_length} ConcordantPair \
		--genomeDir "$rawGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn \$readFilesIn \
		--readFilesCommand gunzip -c \
		--outFilterMultimapNmax 3 \
		--outFileNamePrefix "./" \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--outSAMattrRGline $RG \
		--sjdbGTFfile "$genomeGTF"

	mv ./SJ.out.tab ./${sample}.SJ.out.tab
	mv "./Aligned.out.bam" "./${sample}.pass1.bam"
	mv "./Log.final.out" "./${sample}_Log.final.out"
	"""
}

if(params.umi) {
	// Change read name, the "_" into a ":" before the UMI in read name;
	// Create an unmapped BAM and mapped it with STAR
	// see: https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md
	process umi_stat_and_consensus{

		cpus { params.CPU_umi }
		label 'multicore'
		label 'retriable'
		storeDir { "${params.out}/fgbio" }

		input:
		set val(BAM), val(sample), val(type), val(RG) from BAM_pass1

		output:
		set val(sample), file("${sample}_family_size_histogram.txt") into UMI_stat
		file "${sample}_family_size_histogram.txt" into UMI_table
		set file("${sample}.consensus_R1.fastq.gz"), file("${sample}.consensus_R2.fastq.gz"), val(sample), val(type), val(RG) into FASTQ_STAR2
		set val(sample), file("${sample}.consensus.bam") into BAM_unmapped

		"""
		### Create a temporary directory for fgbio execution
		tmpdir="${sample}_\$(date +%H%M%S)"
		mkdir "\${tmpdir}"

		### fgbio command
		fgBioExe="java -Xmx4g -Djava.io.tmpdir=\${tmpdir} -jar \$fgbio"

		### Function to run at the end to clean the temporary files
		function cleanup()
		{
			rm -rf "${sample}.changeName.bam" "${sample}.copy.bam" "${sample}.sort.bam" "${sample}.mate.bam" "${sample}.grpUmi.bam" "\${tmpdir}"
		}

		### Clean the temporary file when the program exit
		trap cleanup EXIT

		### Change the "_" into a ":" before the UMI in read name
		samtools view -h "${BAM}" | sed -r 's/(^[^\t]*:[0-9]*)_([ATCGN]*)\t/\\1:\\2\t/' | samtools view -b > "${sample}.changeName.bam"

		### Put UMI as a tag in the bam file
		\${fgBioExe} --async-io --compression 1 CopyUmiFromReadName \
			--input="${sample}.changeName.bam" \
			--output="${sample}.copy.bam"

		### Put mate info after sorting
		\${fgBioExe} --async-io --compression 1 SortBam \
			--input="${sample}.copy.bam" \
			--output="${sample}.sort.bam" \
			--sort-order=Queryname
		\${fgBioExe} --async-io --compression 1 SetMateInformation \
			--input="${sample}.sort.bam" \
			--output="${sample}.mate.bam"

		### Group reads per UMI
		\${fgBioExe} --async-io --compression 1 GroupReadsByUmi \
			--input="${sample}.mate.bam" \
			--output="${sample}.grpUmi.bam" \
			--family-size-histogram="${sample}_family_size_histogram.txt" \
			--raw-tag=RX \
			--assign-tag=MI \
			--strategy=Adjacency \
			--edits=1

		### Get consensus
		\${fgBioExe} --async-io CallMolecularConsensusReads \
			--input="${sample}.grpUmi.bam" \
			--output="${sample}.consensus.bam" \
			--error-rate-pre-umi 30 \
			--error-rate-post-umi 30 \
			--output-per-base-tags true \
			--min-reads 1 \
			--max-reads 50 \
			--min-input-base-quality 10 \
			--read-name-prefix="csr" \
			--threads "${params.CPU_umi}"

		### Convert into FASTQ
		samtools collate -u -O "${sample}.consensus.bam" | \
			samtools fastq \
			-1 "${sample}.consensus_R1.fastq.gz" \
			-2 "${sample}.consensus_R2.fastq.gz" \
			-n
		"""
	}
}
else {
	// Bypass
	FASTQ_STAR1_copy.set{ FASTQ_STAR2 }
}

// Get the UMI duplication stat in the FASTQC
if(params.umi) {
	process umi_plot {

		cpus 1
		label 'retriable'
		storeDir { "${params.out}/QC/umi" }

		input:
		set val(sample), file(umiHist) from UMI_stat
		file umi_fastqc from file("${baseDir}/scripts/umi_stat.R")

		output:
		file "${sample}_mqc.yaml" into QC_umi

		"""
		Rscript --vanilla "$umi_fastqc" "$sample" "$umiHist"
		"""
	}
} else {
	QC_umi = Channel.value(file("$baseDir/in/dummy.tsv"))
}

// Get the UMI table
if(params.umi) {
	process umi_table {

		cpus 1
		label 'retriable'
		storeDir { "${params.out}/QC/umi" }

		input:
		file "*" from UMI_table.collect()
		file umi_table from file("${baseDir}/scripts/umi_table.R")

		output:
		file "umi_table_mqc.yaml" into QC_umi_table

		"""
		Rscript --vanilla "$umi_table"
		"""
	}
} else {
	QC_umi_table = Channel.value(file("$baseDir/in/dummy.tsv"))
}

// Build a new genome from STAR pass 1
process STAR_reindex {
	
	cpus 2
	label 'multicore'
	label 'retriable'
	storeDir { params.out }
	
	input:
	file SJ from SJ_pass1.collect()
	file R1 from file("$baseDir/in/dummy_R1.fastq")
	file R2 from file("$baseDir/in/dummy_R2.fastq")
	file genomeGTF from genomeGTF
	file rawGenome from rawGenome_reindex
	
	output:
	file("${params.genome}_${params.title}") into reindexedGenome
	
	"""
	mkdir -p "./reindex"
	STAR \
		--runThreadN 2 \
		--genomeDir "$rawGenome" \
		--readFilesIn "$R1" "$R2" \
		--sjdbFileChrStartEnd $SJ \
		--limitSjdbInsertNsj 5000000 \
		--sjdbInsertSave All \
		--sjdbGTFfile "$genomeGTF" \
		--outFileNamePrefix "./reindex/" \
		--outSAMtype None
	mv "./reindex/Log.out" "./reindex/_STARgenome/"
	mv ./reindex/_STARgenome/ ./${params.genome}_${params.title}/
	"""
}

// STAR second pass
// TODO shared-memory
process STAR_pass2 {
	
	cpus { params.CPU_align2 }
	label 'multicore'
	label 'retriable'
	storeDir { "${params.out}/STAR_pass2" }
	
	input:
	set file(R1), file(R2), val(sample), val(type), val(RG) from FASTQ_STAR2
	file reindexedGenome
	file genomeGTF from genomeGTF

	output:
	set val(sample), val(type), file("${sample}.DNA.temp.bam") into genomic_temp_BAM
	set val(sample), val(type), file("${sample}.RNA.bam") optional true into transcriptomic_BAM
	set val(sample), val(type), file("${sample}.isize.txt") into isize_sample
	file "${sample}_SJ.out.tab" into junctions_STAR
	file "${sample}_Chimeric.out.junction" into chimeric_STAR
	file "${sample}.isize.txt" into isize_table
	file "${sample}_Log.final.out" into QC_STAR_pass2
	
	"""
	# FASTQ files
	if [ "$type" = "paired" ];   then readFilesIn="\\"${R1.join(",")}\\" \\"${R2.join(",")}\\""
	elif [ "$type" = "single" ]; then readFilesIn="\\"${R1.join(",")}\\""
	else                         echo "Unknow type '$type'"; exit 1
	fi
	
	# Align
	mkdir -p "./$sample"
	STAR \
		--runThreadN ${params.CPU_align2} \
		--twopassMode None \
		--genomeDir "$reindexedGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn \$readFilesIn \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "./${sample}/" \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--chimOutType Junctions WithinBAM \
		--quantMode TranscriptomeSAM \
		--outSAMattrRGline $RG \
		--sjdbGTFfile "$genomeGTF" \
		--alignEndsProtrude ${params.umi_length} ConcordantPair \
		--alignInsertionFlush Right \
		--alignSJDBoverhangMin 4 \
		--alignSJstitchMismatchNmax 3 -1 3 3 \
		--alignSplicedMateMapLmin 16 \
		--alignSplicedMateMapLminOverLmate 0 \
		--chimJunctionOverhangMin 8 \
		--chimScoreJunctionNonGTAG -4 \
		--chimSegmentMin 10 \
		--chimMultimapNmax 1 \
		--chimNonchimScoreDropMin 10 \
		--outFilterMultimapNmax 3 \
		--outFilterMismatchNmax 5 \
		--outSJfilterOverhangMin 16 8 8 8 \
		--outSJfilterDistToOtherSJmin 0 0 0 0 \
		--peOverlapNbasesMin 12 \
		--peOverlapMMp 0.1

	mv "./${sample}/Log.final.out" "./${sample}_Log.final.out"
	mv "./${sample}/SJ.out.tab" "./${sample}_SJ.out.tab"
	mv "./${sample}/Chimeric.out.junction" "./${sample}_Chimeric.out.junction"
	mv "./${sample}/Aligned.out.bam" "./${sample}.DNA.temp.bam"
	mv "./${sample}/Aligned.toTranscriptome.out.bam" "./${sample}.RNA.bam"
	
	# Export ISIZE sample (empty in single-end)
	samtools view -f 0x2 -f 0x80 "./${sample}.RNA.bam" | cut -f9 | head -1000000 > "./${sample}.isize.txt"
	"""
}

// Get the median insert size per sample
process insertSize_table {

	cpus 1
	label 'retriable'
	storeDir { "${params.out}/QC/insertSize" }

	input:
	file "*" from isize_table.collect()
	file insert_table from file("${baseDir}/scripts/insert_table.R")

	output:
	file "isize_table_mqc.yaml" into median_isize_table

	"""
	Rscript --vanilla "$insert_table"
	"""
}

// Prepare FASTA satellite files as requested by GATK
process indexFASTA {

	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { params.store }
	scratch { params.scratch }

	input:
	file genomeFASTA from genomeFASTA

	output:
	set file(genomeFASTA), file("${genomeFASTA.getBaseName()}.dict"), file("${genomeFASTA}.fai") into (indexedFASTA_splitN, indexedFASTA_BQSR, indexedFASTA_Mutect2, indexedFASTA_MergeBamAlignment)

	"""
	# Dictionnary
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CreateSequenceDictionary \
		REFERENCE="$genomeFASTA" \
		OUTPUT="${genomeFASTA.getBaseName()}.dict"
	# Index
	samtools faidx "$genomeFASTA"
	"""
}

if(params.umi) {
	// Merge mapped and unmapped BAM and filter
	process merge_filterBam {

		cpus 2
		label 'retriable'
		storeDir { "${params.out}/mergeBam" }

		input:
		set val(sample), val(type), file(BAM_mapped), file(BAM_unmapped), file(BAM_forUnmappedRead) from genomic_temp_BAM.join(BAM_unmapped.join(BAM_forUnmappedRead))
		set file(genomeFASTA), file(genomeFASTAdict), file(genomeFASTAindex) from indexedFASTA_MergeBamAlignment

		output:
		set val(sample), val(type), file("${sample}.DNA.bam") into genomic_BAM

		"""
		### fgbio command
		fgBioExe="java -Xmx4g -jar \$fgbio"

		### Sort the BAM by query name for gatk - VW need to be put before in STAR_pass2
		mkdir tmp
		java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" SortSam \
			-INPUT "${BAM_mapped}" \
			-OUTPUT mapped_sorted.bam \
			-SORT_ORDER queryname \
			--TMP_DIR "\$(pwd)/tmp"
		java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" SortSam \
			-INPUT "${BAM_unmapped}" \
			-OUTPUT unmapped_sorted.bam \
			-SORT_ORDER queryname \
			--TMP_DIR "\$(pwd)/tmp"

		### Values from ROCHE pipeline
		gatk MergeBamAlignment \
			--ATTRIBUTES_TO_RETAIN X0 \
			--ATTRIBUTES_TO_RETAIN RX \
			--ALIGNED_BAM mapped_sorted.bam \
			--UNMAPPED_BAM unmapped_sorted.bam \
			--OUTPUT "${sample}.temp.bam" \
			--REFERENCE_SEQUENCE "${genomeFASTA}" \
			--SORT_ORDER coordinate \
			--ADD_MATE_CIGAR true \
			--MAX_INSERTIONS_OR_DELETIONS -1 \
			--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
			--ALIGNER_PROPER_PAIR_FLAGS true \
			--CLIP_OVERLAPPING_READS false \
			--TMP_DIR "\$(pwd)/tmp"

		### Values from Thomas (HCL_nextflow)
		samtools sort -n -u "${sample}.temp.bam" |
			\${fgBioExe} --async-io --compression 1 FilterConsensusReads \
			--ref=${genomeFASTA} \
			--input=/dev/stdin \
			--output="${sample}.filterConsensus.bam" \
			--min-reads=1 \
			--max-read-error-rate=0.05 \
			--max-base-error-rate=0.1 \
			--min-base-quality=1 \
			--max-no-call-fraction=0.2

		### Get unmapped read from STAR_pass1 and put them after (from: https://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/)
		### 0x4 4  UNMAP        0x108 264 MUNMAP,SECONDARY
		### 0x8 8  MUNMAP       0x104 260 UNMAP,SECONDARY
		### 0xc 12 UNMAP,MUNMAP 0x100 256 SECONDARY
		samtools view -b -f4 -F264  "${BAM_forUnmappedRead}" > tmps1.bam
		samtools view -b -f8 -F260  "${BAM_forUnmappedRead}" > tmps2.bam
		samtools view -b -f12 -F256 "${BAM_forUnmappedRead}" > tmps3.bam
		samtools merge -o "${sample}.DNA.bam" "${sample}.filterConsensus.bam" "tmps"?".bam"
		"""
	}
} else {
	process skip_merge_filterBam {
		cpus 1
		label 'monocore'
		label 'nonRetriable'
		storeDir { "${params.out}/mergeBam" }
		executor 'local'

		input:
		set val(sample), val(type), file(BAM_mapped) from genomic_temp_BAM

		output:
		set val(sample), val(type), file("${sample}.DNA.bam") into genomic_BAM

		"""
		mv "${BAM_mapped}" "${sample}.DNA.bam"
		"""
	}
}

// Picard MarkDuplicates (mark only, filter later)
// FIXME use as many CPUs as available, whatever the options
// FIXME add a short @PG line (default adds to all reads and mess up with samtools other @PG)
process markDuplicates {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/markDuplicates" }
	scratch { params.scratch }
	
	input:
	set val(sample), val(type), file(BAM) from genomic_BAM
	
	output:
	file "${sample}.txt" into QC_markDuplicates
	set val(sample), val(type), file("${BAM.getBaseName()}.MD.bam") into BAM_marked
	file "${BAM.getBaseName()}.MD.clean" into markDuplicates_clean
	
	"""
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" MarkDuplicates \
		TMP_DIR="." \
		INPUT="$BAM" \
		OUTPUT="${BAM.getBaseName()}.MD.bam" \
		METRICS_FILE="${sample}.txt" \
		ASSUME_SORT_ORDER="queryname" \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_SEQUENCING_DUPLICATES="false" \
		REMOVE_DUPLICATES="false" \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=50 \
		PROGRAM_RECORD_ID=null
	
	# Link input BAM for cleaning
	ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.MD.clean"
	"""
}

// Genomically sort and index
process BAM_sort {
	
	cpus 4
	label 'multicore'
	label 'retriable'
	storeDir { "${params.out}/BAM" }
	
	input:
	set val(sample), val(type), file(BAM) from BAM_marked
	
	output:
	set val(sample), val(type), file("${BAM.getBaseName()}.sort.bam"), file("${BAM.getBaseName()}.sort.bai") into (BAM_sorted, BAM_rnaSeqMetrics, BAM_featureCounts, BAM_secondary, BAM_softClipping)
	file "${BAM.getBaseName()}.sort.bam" into (BAM_splicing, BAM_dup2)
	file "${BAM.getBaseName()}.sort.bai" into BAI_splicing
	file "${BAM.getBaseName()}.sort.clean" into BAM_sort_clean

	"""
	# Sort
	samtools sort -o "${BAM.getBaseName()}.sort.bam" -T ./${sample} -@ 3 "$BAM"
	
	# Index
	samtools index "${BAM.getBaseName()}.sort.bam"
	mv "${BAM.getBaseName()}.sort.bam.bai" "${BAM.getBaseName()}.sort.bai"
	
	# Link input BAM for cleaning
	ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.sort.clean"
	"""
}

if(params.umi) {
	process duplication_UMI_based {

		cpus 1
		label 'monocore'
		label 'retriable'
		storeDir { "${params.out}/QC" }

		input:
		file '*' from BAM_dup1.collect()
		file '*' from BAM_dup2.collect()
		file run_dup from file("${baseDir}/scripts/duplication_umi.sh")

		output:
		file 'duplication_umi.yaml' into dup_umi

		"""
		bash "$run_dup" "${params.out}"
		"""
	}
} else {
	dup_umi = Channel.value(file("$baseDir/in/dummy.tsv"))
}

// Filter out duplicated read, based on a previous MarkDuplicates run
process filterDuplicates {

	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/filterDuplicates" }
	scratch { params.scratch }
	
	when:
	params.varcall
	
	input:
	set val(sample), val(type), file(BAM), file(BAI) from BAM_sorted
	
	output:
	set val(sample), val(type), file("${BAM.getBaseName()}.filter.bam"), file("${BAM.getBaseName()}.filter.bai") into BAM_filtered
	
	"""
	# Filter
	samtools view -b -F 0x400 "$BAM" -o "${BAM.getBaseName()}.filter.bam"
	
	# Index
	samtools index "${BAM.getBaseName()}.filter.bam"
	mv "${BAM.getBaseName()}.filter.bam.bai" "${BAM.getBaseName()}.filter.bai"
	"""
}

/*
// Prepare a raw vcf.gz file downloaded from COSMIC
process COSMIC {

	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { params.store }
	scratch { params.scratch }
	
	when:
	params.varcall
	
	input:
	file COSMIC from rawCOSMIC
	
	output:
	set file("${COSMIC.getBaseName()}.bgz"), file("${COSMIC.getBaseName()}.bgz.tbi") into COSMIC
	
	// FIXME : edit the "contig" list in header as well or BQSR will fail
	
	"""
	# Add 'chr' to chromosome names and bgzip
	gunzip --stdout "$COSMIC" | awk '/^#/ { print \$0; next } { print "chr"\$0 }' | bgzip --stdout > "${COSMIC.getBaseName()}.bgz"
	
	# Create index
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" IndexFeatureFile -I "${COSMIC.getBaseName()}.bgz"
	"""
}
*/

// Picard SplitNCigarReads (split reads with intron gaps into separate reads)
process splitN {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/SplitNCigarReads" }
	scratch { params.scratch }
	
	input:
	set file(genomeFASTA), file(genomeFASTA_dict), file(genomeFASTA_fai) from indexedFASTA_splitN
	set val(sample), val(type), file(BAM), file(BAI) from BAM_filtered
	
	output:
	set val(sample), val(type), file("${BAM.getBaseName()}.splitN.bam"), file("${BAM.getBaseName()}.splitN.bai") into BAM_splitN
	file "${BAM.getBaseName()}.splitN.clean" into splitN_clean
	
	"""
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" SplitNCigarReads \
		--input "$BAM" \
		--reference "$genomeFASTA" \
		--output "${BAM.getBaseName()}.splitN.bam" \
		--tmp-dir "."
	
	# Link input BAM for cleaning
	ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.splitN.clean"
	"""
}

// Compute and apply GATK Base Quality Score Recalibration model
process BQSR {

	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/BQSR" }
	scratch { params.scratch }
	
	input:
	set file(genomeFASTA), file(genomeFASTA_dict), file(genomeFASTA_fai) from indexedFASTA_BQSR
	set file(gnomAD), file(gnomAD_index) from gnomAD_BQSR
	set file(COSMIC), file(COSMIC_index) from COSMIC
	set val(sample), val(type), file(BAM), file(BAI) from BAM_splitN
	
	output:
	set val(sample), val(type), file("${BAM.getBaseName()}.BQSR.bam"), file("${BAM.getBaseName()}.BQSR.bai") into BAM_BQSR
	file "${BAM.getBaseName()}.BQSR.clean" into BQSR_clean
	
	"""
	# Compute model
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" BaseRecalibrator \
		--input "$BAM" \
		--reference "$genomeFASTA" \
		--known-sites "$gnomAD" \
		--known-sites "$COSMIC" \
		--output "${sample}.BQSR" \
		--tmp-dir "."
	
	# Apply model
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" ApplyBQSR \
		--input "$BAM" \
		--reference "$genomeFASTA" \
		--bqsr-recal-file "${sample}.BQSR" \
		--output "${BAM.getBaseName()}.BQSR.bam" \
		--tmp-dir "."
	
	# Link input BAM for cleaning
	ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.BQSR.clean"
	"""
}

// Call variants with GATK Mutect2 (FIXME : --panel-of-normals pon.vcf.gz)
process Mutect2 {

	cpus { params.CPU_mutect }
	label 'multicore'
	label 'nonRetriable'
	storeDir { "${params.out}/Mutect2" }
	scratch { params.scratch }
	
	input:
	set file(genomeFASTA), file(genomeFASTA_dict), file(genomeFASTA_fai) from indexedFASTA_Mutect2
	set file(gnomAD), file(gnomAD_index) from gnomAD_Mutect2
	set val(sample), val(type), file(BAM), file(BAI) from BAM_BQSR
	
	output:
	set file("${sample}.filtered.vcf.gz"), file("${sample}.filtered.vcf.gz.tbi") into filtered_VCF
	set file("${sample}.unfiltered.vcf.gz"), file("${sample}.unfiltered.vcf.gz.tbi") into unfiltered_VCF
	file("${sample}.unfiltered.vcf.gz.stats") into Mutect2_stats
	
	"""
	# Genomic subset
	if [ "${params.window}" = "" ]; then interval=""
	else                                 interval="--intervals ${params.window}"
	fi
	
	# Call variants
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" Mutect2 \$interval \
		--input "$BAM" \
		--reference "$genomeFASTA" \
		--output "${sample}.unfiltered.vcf.gz" \
		--germline-resource "$gnomAD" \
		--native-pair-hmm-threads ${params.CPU_mutect} \
		--independent-mates
	
	# Filter variants
	gatk FilterMutectCalls \
		--variant "${sample}.unfiltered.vcf.gz" \
		--reference "$genomeFASTA" \
		--output "${sample}.filtered.vcf.gz"	
	"""
}

// Prepare refFlat file for Picard
process refFlat {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { params.store }
	
	input:
	file targetGTF from targetGTF
	file gtfToRefFlat from file("${baseDir}/scripts/gtfToRefFlat.R")
	
	output:
	file "${targetGTF}.refFlat" into target_refFlat
	
	"""
	Rscript --vanilla "$gtfToRefFlat" "$targetGTF" "${targetGTF}.refFlat"
	"""	
}

// Prepare rRNA interval list file for Picard
process rRNA_interval {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { params.store }
	
	input:
	file targetGTF from targetGTF
	file chrNameLength from rawGenome_chrom
	
	output:
	file("${targetGTF}.rRNA") into target_rRNA
	
	"""
	# Header (consider unsorted to be safe)
	echo -e "@HD\tVN:1.0\tSO:unsorted" > "./${targetGTF}.rRNA"

	# Chromosomes (from STAR genome)
	sed -r 's/^(.+)\t(.+)\$/@SQ\tSN:\\1\tLN:\\2/' "${chrNameLength}" >> "./${targetGTF}.rRNA"

	# BED-like content
	grep -E 'transcript_(bio)?type "rRNA"' "$targetGTF" | awk -F "\t" '\$3 == "transcript" { id=gensub(/^.+transcript_id \"([^\"]+)\";.+\$/, "\\\\1", "g", \$9); print \$1"\t"\$4"\t"\$5"\t"\$7"\t"id }' >> "./${targetGTF}.rRNA"
	"""
}

// Picard's CollectRnaSeqMetrics
process rnaSeqMetrics {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { "${params.out}/QC/rnaSeqMetrics" }

	input:
	set val(sample), val(type), file(BAM), file(BAI) from BAM_rnaSeqMetrics
	file target_rRNA
	file target_refFlat
	
	output:
	file "${sample}.RNA_Metrics" into QC_rnaSeqMetrics
	
	"""
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CollectRnaSeqMetrics \
		INPUT=$BAM \
		OUTPUT="./${sample}.RNA_Metrics" \
		REF_FLAT="$target_refFlat" \
		RIBOSOMAL_INTERVALS="$target_rRNA" \
		STRAND_SPECIFICITY="${params.stranded_Picard}" \
		ASSUME_SORTED=true
	"""
}

// Count reads in transcripts using featureCounts
process featureCounts {
	
	cpus 2
	label 'multicore'
	label 'retriable'
	storeDir { "${params.out}/featureCounts" }
	
	input:
	set val(sample), val(type), file(BAM), file(BAI) from BAM_featureCounts
	file targetGTF from targetGTF
	
	output:
	file "annotation.tsv" into featureCounts_annotation
	file "${sample}_counts.rds" into featureCounts_counts
	file "${sample}_stats.rds" into featureCounts_stats
	
	"""
	#!/usr/bin/env Rscript --vanilla
	
	# Dependency
	library(Rsubread)
	
	# Count reads in genes (~5 minutes with 8 threads)
	dir.create("./tmp")
	out <- featureCounts(
		files = "$BAM",
		annot.ext = "$targetGTF",
		isGTFAnnotationFile = TRUE,
		allowMultiOverlap = FALSE,
		minMQS = 0L,
		strandSpecific = ${params.stranded_Rsubread},
		isPairedEnd = ("$type" == "paired"),
		requireBothEndsMapped = TRUE,
		autosort = FALSE,
		nthreads = 2,
		tmpDir = "./tmp"
	)
	file.remove("./tmp")
	
	# Export annotation (once would be enough)
	write.table(out\$annotation, file="./annotation.tsv", row.names=FALSE, quote=FALSE, sep="\t")
	
	# Export counts
	counts <- out\$counts
	colnames(counts) <- "${sample}"
	saveRDS(counts, file="./${sample}_counts.rds")
	
	# Export counts
	stats <- out\$stat
	colnames(stats)[2] <- "${sample}"
	saveRDS(stats, file="./${sample}_stats.rds")
	"""
}

// Use edgeR to compute QC
process edgeR {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/edgeR" }
	
	when:
	params.finalize
	
	input:
	file 'annotation' from featureCounts_annotation.first()
	file 'countFiles' from featureCounts_counts.collect()
	file edgeR from file("${baseDir}/scripts/edgeR.R")
	
	output:
	file 'counts.tsv' into counts
	file 'CPM.tsv' into CPM
	file 'RPK.tsv' into RPK
	file 'edgeR.yaml' into QC_edgeR_general
	file 'edgeR_mqc.yaml' into QC_edgeR_section
	
	"""
	Rscript --vanilla "$edgeR" "$annotation" "." $countFiles
	"""	
}

// Estimate insert size distribution
process insertSize {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { "${params.out}/QC/insertSize" }
	
	when:
	type == "paired"
	
	input:
	set val(sample), val(type), file(isize) from isize_sample
	file insertSize from file("${baseDir}/scripts/insertSize.R")
	
	output:
	file "${sample}_mqc.yaml" into QC_insert
	
	"""
	Rscript --vanilla "$insertSize" "$sample" "$isize" > "./${sample}_mqc.yaml"
	"""
}

// Quantify secondary alignments with SAMtools
// TODO : general stats
process secondary {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { "${params.out}/QC/secondary" }
	
	input:
	set val(sample), val(type), file(BAM), file(BAI) from BAM_secondary
	
	output:
	file "${sample}_mqc.yaml" into QC_secondary
	
	"""
	# Total read count (fast, from index)
	total=\$(samtools idxstats $BAM | awk 'BEGIN{x=0} {x=x+\$3+\$4} END{print x}')
	
	# Unaligned
	una=\$(samtools view -c -f 0x4 $BAM)
	
	# Secondary alignments
	sec=\$(samtools view -c -f 0x100 $BAM)
	
	# Primary alignments (deduced)
	pri=\$(bc <<< "scale=2; \$total-\$una-\$sec")
	
	# MultiQC regular file header
	cat <<-EOF > "./${sample}_mqc.yaml"
		id: 'Secondary_section'
		section_name: 'Secondary alignments'
		description: 'as a proportion of all alignments (reads) returned by STAR. A read aligning in multiple locations is duplicated in the BAM file (one entry for each alignment), the proportion of secondary alignments shows to which extent the read count was artificially increased by this phenomenon.'
		plot_type: 'bargraph'
		pconfig:
		    id: 'Secondary_bargraph'
		    title: 'Secondary alignments'
		data:
		    ${sample}: {Primary: \${pri}, Unaligned: \${una}, Secondary: \${sec}}
	EOF
	"""
}

// Plot soft-clipping lengths on read ends
// TODO : general stats
process softClipping {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { "${params.out}/QC/softClipping" }
	
	input:
	set val(sample), val(type), file(BAM), file(BAI) from BAM_softClipping
	file softClipping from file("${baseDir}/scripts/softClipping.R")
	
	output:
	file "${sample}_*_mqc.yaml" into QC_softClipping
	
	"""
	Rscript --vanilla "$softClipping" "$sample" "$BAM"
	"""
}

// Collect QC files into a single report
process MultiQC {

	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/QC" }

	when:
	params.finalize

	input:
	file conf from file("$baseDir/in/multiqc.conf")
	file 'edgeR.yaml' from QC_edgeR_general
	file 'edgeR_mqc.yaml' from QC_edgeR_section
	file 'STAR_pass1/*' from QC_STAR_pass1.collect()
	file 'STAR_pass2/*' from QC_STAR_pass2.collect()
	file 'FastQC/raw/*_fastqc.zip' from QC_FASTQC_raw.collect()
	file 'FastQC/trimmed/*_fastqc.zip' from QC_FASTQC_trimmed.collect()
	file 'markDuplicates/*' from QC_markDuplicates.collect()
	file 'rnaSeqMetrics/*' from QC_rnaSeqMetrics.collect()
	file 'insertSize/*' from insertSize_bypass.mix(QC_insert).collect()
	file 'secondary/*' from QC_secondary.collect()
	file 'softClipping/*' from QC_softClipping.collect()
	file 'umi/*_mqc.yaml' from QC_umi.collect()
	file 'umi_table_mqc.yaml' from QC_umi_table
	file 'isize_table_mqc.yaml' from median_isize_table
	file 'cutadapt/*' from QC_cutadapt.collect()
	file 'duplication_umi.yaml' from dup_umi

	output:
	file "${params.MQC_title}_multiqc_report_data.zip" into MultiQC_data
	file "${params.MQC_title}_multiqc_report.html" into MultiQC_report

	"""
	multiqc --title "${params.MQC_title}" --comment "${params.MQC_comment}" --outdir "." --config "${conf}" --config "./edgeR.yaml" --config "./umi_table_mqc.yaml" --config "./duplication_umi.yaml" --config "./isize_table_mqc.yaml" --zip-data-dir --interactive --force "."
	"""
}

// Remove unnecessary BAM (unstorable process)
process clean_DNA_BAM {
	
	cpus 1
	label 'nonRetriable'
	
	// Never scratch
	scratch false
	stageInMode 'symlink'
	executor 'local'
	
	when:
	params.clean_BAM
	
	input:
	file(clean) from BAM_sort_clean.mix(markDuplicates_clean, splitN_clean, BQSR_clean)
	
	"""
	# BAM file
	BAM="\$(readlink "$clean")"
	echo -n '' > \$BAM
	
	# BAI file
	BAI="\${BAM%.bam}.bai"
	if [ -f "\$BAI" ]; then echo -n '' > \$BAI; fi
	"""
}

// Remove unnecessary BAM (unstorable process)
process clean_RNA_BAM {
	
	cpus 1
	label 'nonRetriable'
	
	// Never scratch
	scratch false
	stageInMode 'symlink'
	executor 'local'
	
	when:
	params.clean_BAM && !params.RNA_BAM
	
	input:
	set val(sample), val(type), file(BAM) from transcriptomic_BAM
	
	"""
	echo -n '' > "\$(readlink "$BAM")"
	"""
}

// Prepare introns and exon track files
process annotation {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { params.store }
	
	when:
	params.splicing
	
	input:
	file genomeGTF from genomeGTF
	file script from file("${baseDir}/scripts/annotation.R")
	
	output:
	file "${genomeGTF}.introns.rds" into introns
	file "${genomeGTF}.exons.rdt" into (exons_collect, exons_filter)
	
	"""
	Rscript --vanilla "$script" "$genomeGTF" "$params.species" "$params.genome" "$params.chromosomes"
	"""
}

// Collect all splicing events
process splicing_collect {
	
	cpus { params.CPU_splicing }
	label 'multicore'
	label 'nonRetriable'
	storeDir { "${params.out}/splicing" }
	
	when:
	params.splicing
	
	input:
	file exons from exons_collect
	file introns from introns
	file 'junctionFiles/*' from junctions_STAR.collect()
	file 'chimericFiles/*' from chimeric_STAR.collect()
	file script from file("${baseDir}/scripts/splicing_collect.R")
	file "transcripts.tsv" from transcripts
	
	output:
	file("*.rds") into splicing_events
	
	"""
	Rscript --vanilla "$script" ${params.CPU_splicing} "$exons" "$introns" "$params.chromosomes" $params.min_reads_unknown "transcripts.tsv"
	"""
}

// Collect all splicing events
process splicing_filter {
	
	cpus 1
	label 'monocore'
	label 'nonRetriable'
	storeDir { "${params.out}/splicing" }
	
	when:
	params.splicing
	
	input:
	file exons from exons_filter
	file '*' from splicing_events
	file script from file("${baseDir}/scripts/splicing_filter.R")
	file '*' from BAM_splicing.collect()
	file '*' from BAI_splicing.collect()
	file targetGTF from targetGTF
	
	output:
	file("I-${params.min_I}_PSI-${params.min_PSI}_${params.symbols.take(50)}(${params.symbols.split(',').size()})_${params.classes}_${params.focus.replaceAll(':','-')}") into splicing_output
	file("depth") into splicing_depth
	
	"""
	Rscript --vanilla "$script" ${params.CPU_splicing} "$targetGTF" "$exons" ${params.plot} ${params.min_I} ${params.min_PSI} "$params.symbols" "$params.classes" "$params.focus"
	"""
}
