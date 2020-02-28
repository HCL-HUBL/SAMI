#!/usr/bin/env nextflow

/*
 * RNA-seq pipeline
 * <sylvain.mareschal@lysarc.org>
 *
 * nextflow run RNA-seq.nf -with-singularity /dev/shm/RNA-seq.sif --title 'Test' --FASTQ 'data/test' --stranded 'R2' --RG_CN 'Integragen' --RG_PL 'ILLUMINA' --RG_PM 'HiSeq2000' --CPU_index 48 --CPU_align1 6 --CPU_align2 16
 */

// Run characteristics (no default value)
params.FASTQ = ''
params.stranded = ''
params.RG_CN = ''
params.RG_PL = ''
params.RG_PM = ''
params.title = ''

// CPU to use (no default value)
params.CPU_index = 0
params.CPU_align1 = 0
params.CPU_align2 = 0

// Mandatory values
if(params.FASTQ == '')                      error "ERROR: --FASTQ must be provided"
if(params.CPU_index <= 0)                   error "ERROR: --CPU_index must be a positive integer (suggested: all available CPUs)"
if(params.CPU_align1 <= 0)                  error "ERROR: --CPU_align1 must be a positive integer (suggested: 6+)"
if(params.CPU_align2 <= 0)                  error "ERROR: --CPU_align2 must be a positive integer (suggested: 6+)"
if(params.title == '')                      error "ERROR: --title must be provided"
if(params.title ==~ /.*[^A-Za-z0-9_\.-].*/) error "ERROR: --title can only contain letters, digits, '.', '_' or '-'"

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

// STAR index (files not provided)
params.species = 'Human'
params.genome = 'GRCh38'
params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
params.genomeFASTA = ''   /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz */
params.genomeGTF = ''     /* ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz */

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
params.maxRetries = 2

// Whether to handle single-end data (R1 only) or consider missing R2 file as an error
params.single = false



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
	sampleDirectory.eachFileMatch(~/(?i).*1.f(ast)?q.gz/) { R1_file ->
		// Corresponding R2 file
		R2_name = R1_file.name.replaceFirst(/(?i)1(.f(ast)?q\.gz)$/, '2$1')
		R2_file = file("${params.FASTQ}/${sample}/${R2_name}")
		if(R2_file.exists()) {
			// Use R2 as R2
			anyPE = true
		} else {
			// Corresponding R3 file
			R3_name = R1_file.name.replaceFirst(/(?i)1(.f(ast)?q\.gz)$/, '3$1')
			R3_file = file("${params.FASTQ}/${sample}/${R3_name}")
			
			if(R3_file.exists()) {
				// Use R3 as R2
				R2_file = R3_file;
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
genomeGTF = Channel.value(file(params.genomeGTF))
genomeFASTA = Channel.value(file(params.genomeFASTA))
headerRegex = Channel.value(file("$baseDir/in/FASTQ_headers.txt"))



// Build RG line from 1st read of each FASTQ file pair bundle
process FASTQ {
	
	cpus 1
	
	// Never scratch to avoid full copy of output in ram-disk
	scratch false
	stageInMode 'symlink'
	executor 'local'
	
	input:
	set file(R1), file(R2), val(sample), val(type) from FASTQ
	file regex from headerRegex
	
	output:
	set file(R1), file(R2), val(sample), val(type), stdout into (FASTQ_STAR1, FASTQ_STAR2)
	file(R1) into FASTQ_R1
	file(R2) into FASTQ_R2
	
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
		
		# RG definition for BAM header
		RG <- c(
			RG,
			sprintf(
				"ID:%s_%i%s\tCN:%s\tPL:%s\tPM:%s\tPU:%s\tSM:%s",
				"${sample}",
				i,
				BC,
				"${params.RG_CN}",
				"${params.RG_PL}",
				"${params.RG_PM}",
				paste(ID1, collapse=":"),
				"${sample}"
			)
		)
	}
	
	# Print final RG to stdout
	cat(paste(RG, collapse=" , "))
	"""
}

// Run FastQC on individual FASTQ files
process FastQC {
	
	cpus 1
	label 'monocore'
	label 'retriable'

	storeDir { "${params.out}/QC/FastQC" }
	
	when:
	!(FASTQ.name ==~ /input\.[0-9]+/)
	
	input:
	file FASTQ from FASTQ_R1.mix(FASTQ_R2).flatten()
	file adapters from file("$baseDir/in/adapters.tab")
	
	output:
	file "${FASTQ.name.replaceFirst(/\.gz$/, '').replaceFirst(/\.f(ast)?q$/, '') + '_fastqc.zip'}" into QC_FASTQC
	
	"""
	fastqc $FASTQ --adapters "$adapters" -o "."
	"""
}

// Build STAR index
// 2019-08-28 CALYM : 27% of 48 CPU usage, 40 GB RAM (scratch = false), 44 min
process STAR_index {
	
	cpus { params.CPU_index }
	label 'multicore'
	storeDir { params.store }
	
	input:
	file genomeFASTA from genomeFASTA
	file genomeGTF from genomeGTF
	
	output:
	file("${params.genome}_raw") into (rawGenome_pass1, rawGenome_reindex, rawGenome_interval)
	
	"""
	mkdir -p "./${params.genome}_raw"
	STAR \
		--runThreadN ${params.CPU_index} \
		--runMode genomeGenerate \
		--genomeDir "./${params.genome}_raw" \
		--genomeFastaFiles "$genomeFASTA" \
		--sjdbGTFfile "$genomeGTF"
	mv "Log.out" "./${params.genome}_raw"
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
	
	output:
	file("${sample}.SJ.out.tab") into SJ_pass1
	
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
		--genomeDir "$rawGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn \$readFilesIn \
		--readFilesCommand gunzip -c \
		--outFilterMultimapNmax 3 \
		--outFileNamePrefix "./" \
		--outSAMtype None
	mv ./SJ.out.tab ./${sample}.SJ.out.tab
	"""
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
	
	output:
	set val(sample), val(type), file("${sample}.DNA.bam") into genomic_BAM
	set val(sample), val(type), file("${sample}.RNA.bam") optional true into transcriptomic_BAM
	set val(sample), val(type), file("${sample}_SJ.out.tab") into junctions_STAR
	set val(sample), val(type), file("${sample}.isize.txt") into isize_sample
	file "${sample}_Log.final.out" into QC_STAR
	
	"""
	# Abort on error (to avoid cleaning BAM too early)
	set -e
	
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
		--outFilterMultimapNmax 3 \
		--outFileNamePrefix "./${sample}/" \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--outSAMattrRGline $RG \
		--quantMode TranscriptomeSAM
	mv "./${sample}/Log.final.out" "./${sample}_Log.final.out"
	mv "./${sample}/SJ.out.tab" "./${sample}_SJ.out.tab"
	mv "./${sample}/Aligned.out.bam" "./${sample}.DNA.bam"
	mv "./${sample}/Aligned.toTranscriptome.out.bam" "./${sample}.RNA.bam"
	
	# Export ISIZE sample (empty in single-end)
	samtools view -f 0x2 -f 0x80 "./${sample}.RNA.bam" | cut -f9 | head -1000000 > "./${sample}.isize.txt"
	
	# Discard RNA BAM (save disk space in work) if not exported
	if [ "${params.RNA_BAM}" == "false" ] && [ "${params.clean_BAM}" == "true" ]; then rm "./${sample}.RNA.bam"; fi
	"""
	// --chimSegmentMin ...
	// --chimOutType WithinBAM
}

// Genomically sort and index
process BAM_sort {
	
	cpus 4
	label 'multicore'
	label 'retriable'
	storeDir { "${params.out}/BAM" }
	
	input:
	set val(sample), val(type), file(BAM) from genomic_BAM
	
	output:
	set val(sample), val(type), file("${sample}.DNA.sorted.bam"), file("${sample}.DNA.sorted.bam.bai") into BAM_rnaSeqMetrics, BAM_featureCounts, BAM_secondary
	
	"""
	# Abort on error (to avoid cleaning BAM too early)
	set -e
	
	# Get file name from Nextflow
	BAM="$BAM"
	
	# Sort
	samtools sort -o \${BAM%.bam}.sort.bam -T ./${sample} -@ 3 \$BAM
	
	# Index
	samtools index \${BAM%.bam}.sort.bam
	
	# Empty the original unsorted BAM file, following the symlink used for stage-in
	if [ "${params.clean_BAM}" == "true" ]; then echo -n '' > "\$(readlink "$BAM")"; fi
	"""
}

// Prepare refFlat file for Picard
// 2019-08-28 CALYM : 100% of 1 CPU usage, 378 MB RAM (scratch = false), 2.8 min
process refFlat {
	
	cpus 1
	label 'monocore'
	storeDir { params.store }
	
	input:
	file genomeGTF from genomeGTF
	file gtfToRefFlat from file("${baseDir}/scripts/gtfToRefFlat.R")
	
	output:
	file "${genomeGTF}.refFlat" into genomeRefFlat
	
	"""
	Rscript --vanilla "$gtfToRefFlat" "$genomeGTF" "${genomeGTF}.refFlat"
	"""	
}

// Prepare rRNA interval list file for Picard
// 2019-08-28 CALYM : 96% of 1 CPU usage, 23.65 MB RAM (scratch = false), 0 min
process rRNA_interval {
	
	cpus 1
	label 'monocore'
	storeDir { params.store }
	
	input:
	file genomeGTF from genomeGTF
	file rawGenome from rawGenome_interval
	
	output:
	file 'rRNA.interval_list' into rRNA_interval
	
	"""
	# Header (consider unsorted to be safe)
	echo -e "@HD\tVN:1.0\tSO:unsorted" > "./rRNA.interval_list"

	# Chromosomes (from STAR genome)
	sed -r 's/^(.+)\t(.+)\$/@SQ\tSN:\\1\tLN:\\2/' "${rawGenome}/chrNameLength.txt" >> "./rRNA.interval_list"

	# BED-like content
	grep 'transcript_type "rRNA"' "$genomeGTF" | awk -F "\t" '\$3 == "transcript" { id=gensub(/^.+transcript_id \"([^\"]+)\";.+\$/, "\\\\1", "g", \$9); print \$1"\t"\$4"\t"\$5"\t"\$7"\t"id }' >> "./rRNA.interval_list"
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
	file rRNA_interval
	file genomeRefFlat
	
	output:
	file "${sample}.RNA_Metrics" into QC_rnaSeqMetrics
	
	"""
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$Picard" CollectRnaSeqMetrics \
		INPUT=$BAM \
		OUTPUT="./${sample}.RNA_Metrics" \
		REF_FLAT="$genomeRefFlat" \
		RIBOSOMAL_INTERVALS="$rRNA_interval" \
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
	file genomeGTF from genomeGTF
	file genomeRefFlat
	
	output:
	file "annotation.csv" into featureCounts_annotation
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
		annot.ext = "$genomeGTF",
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
	write.csv(out\$annotation, file="./annotation.csv", row.names=FALSE, quote=FALSE)
	
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
	storeDir { "${params.out}/edgeR" }
	
	when:
	params.finalize
	
	input:
	file 'annotation' from featureCounts_annotation.first()
	file 'countFiles' from featureCounts_counts.collect()
	file edgeR from file("${baseDir}/scripts/edgeR.R")
	
	output:
	file 'all_counts.rds' into countMatrix
	file 'edgeR.yaml' into QC_edgeR_general
	file 'edgeR_mqc.yaml' into QC_edgeR_section
	
	"""
	Rscript --vanilla "$edgeR" "$annotation" "." $countFiles
	"""	
}

// Estimate insert size distribution
process insertSize {
	
	cpus 2
	label 'multicore'
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
	
	# Secondary alignments
	sec=\$(samtools view -c -f 0x100 $BAM)
	
	# Primary alignments (deduced)
	pri=\$(bc <<< "scale=2; \$total-\$sec")
	
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
		    ${sample}: {Primary: \${pri}, Secondary: \${sec}}
	EOF
	"""
}

// Collect QC files into a single report
process MultiQC {
	
	cpus 1
	label 'monocore'
	storeDir { "${params.out}/QC" }
	
	when:
	params.finalize
	
	input:
	file conf from file("$baseDir/in/multiqc.conf")
	file 'edgeR.yaml' from QC_edgeR_general
	file 'edgeR_mqc.yaml' from QC_edgeR_section
	file 'STAR/*' from QC_STAR.collect()
	file 'FASTQC/*' from QC_FASTQC.collect()
	file 'rnaSeqMetrics/*' from QC_rnaSeqMetrics.collect()
	file 'insertSize/*' from insertSize_bypass.mix(QC_insert).collect()
	file 'secondary/*' from QC_secondary.collect()
	
	output:
	file "${params.MQC_title}_multiqc_report_data.zip" into MultiQC_data
	file "${params.MQC_title}_multiqc_report.html" into MultiQC_report
	
	"""
	multiqc --title "${params.MQC_title}" --comment "${params.MQC_comment}" --outdir "." --config "${conf}" --config "./edgeR.yaml" --zip-data-dir --interactive --force "."
	"""
}

// Reshape STAR junction file into a Rgb table
process junctions {
	
	cpus 1
	label 'monocore'
	label 'retriable'
	storeDir { "${params.out}/junctions" }
	
	input:
	set val(sample), val(type), file(SJ_tab) from junctions_STAR
	
	output:
	file "${sample}.rdt" into junctions_Rgb
	
	"""
	#!/usr/bin/env Rscript --vanilla
	
	# Dependency
	library(Rgb)
	
	# Parse STAR junction file
	tab <- read.table(
		"${SJ_tab}", sep="\t", quote=NULL, comment.char="",
		col.names = c("chrom", "start", "end", "strand", "motif", "annotated", "reads.uni", "reads.multi", "overhang"),
		colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer"),
	)
	
	# Reshape chromosome
	tab\$chrom <- factor(sub("^chr", "", tab\$chrom), levels=strsplit("${params.chromosomes}", split=",", fixed=TRUE)[[1]])
	tab <- tab[ !is.na(tab\$chrom) ,]

	# Reshape strand
	tab[ tab\$strand == 1L , "strand" ] <- "+"
	tab[ tab\$strand == 2L , "strand" ] <- "-"
	tab[ tab\$strand == 0L , "strand" ] <- NA
	tab\$strand <- factor(tab\$strand, levels=c("-","+"))
	
	# Add read counts
	tab\$reads <- tab\$reads.uni + tab\$reads.multi
	tab <- tab[ tab\$reads >= 3L ,]
	
	# Simplify
	for(k in c("motif", "annotated", "reads.uni", "reads.multi", "overhang")) tab[[k]] <- NULL
	
	# Convert to Rgb track
	tab\$name <- as.character(tab\$reads)
	jun <- track.table(tab, .name="${sample}", .organism="${params.species}", .assembly="${params.genome}")
	jun\$setParam("fillColor", function() { x <- slice\$reads / max(slice\$reads); rgb(red=1-x, green=1-x, blue=1) })
	jun\$setParam("maxElements", 100L)
	
	# Export
	saveRDT(jun, file="./${sample}.rdt")
	"""
}

