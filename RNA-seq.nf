#!/usr/bin/env nextflow

/*
 * RNA-seq pipeline
 * <sylvain.mareschal@lysarc.org> 2019-09-04
 */

/*
// Make test data
for sample in Sample_700650 Sample_700652 Sample_700656
do
   for file in $(ls "/home/sylvain/PRIMA_data/RNAseq/FASTQ/$sample")
   do
      echo $file
      gunzip --stdout "/home/sylvain/PRIMA_data/RNAseq/FASTQ/$sample/$file" | head -100000 | gzip --stdout > "$sample/$file"
   done
done
*/


// FASTQ root
params.FASTQ = 'PRIMA_data'
params.readLength = 76
params.strandness = 'SECOND_READ_TRANSCRIPTION_STRAND'
params.strandSpecificity = '2L'
params.RG_CN = 'Integragen'
params.RG_PL = 'ILLUMINA'
params.RG_PM = 'HiSeq2000'



// Long-term storage
storeDir = "${baseDir}/store"
publishDir = "${baseDir}/out"

// How to deal with output files (link from ./work or move from /dev/shm)
params.debug = true
if(params.debug) {
	scratchMode = 'false'
	publishMode = 'symlink'
} else {
	scratchMode = 'ram-disk'
	publishMode = 'move'
}

// CPU to use (in total and per sample)
params.CPU_STAR_index = 48
params.CPU_STAR_align = 6

// STAR index
params.species = 'Human'
params.genome = 'GRCh38'
params.chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
params.genomeFASTA = '/data/P71732_DATA01_SAS/tools/local_install/STAR_GRCh38/GRCh38.primary_assembly.genome.fa'
params.genomeGTF = '/data/P71732_DATA01_SAS/tools/local_install/STAR_GRCh38/gencode.v29.annotation.gtf'

// Tools
params.STAR = '/data/P71732_DATA01_SAS/tools/local_install/STAR-2.6.1d/bin/Linux_x86_64/STAR'
params.samtools = '/data/P71732_DATA01_SAS/tools/local_install/samtools-1.9/samtools'
params.fastqc = '/data/P71732_DATA01_SAS/tools/local_install/FastQC-0.11.8/fastqc'
params.picard = '/data/P71732_DATA01_SAS/tools/local_install/picard-2.18.20/picard.jar'
params.R_packages = '/data/P71732_DATA01_SAS/tools/local_install/R-packages'

// Multi-QC annotation
params.MQC_title = 'RNA-seq analysis'
params.MQC_comment = 'Processed through RNA-seq.nf'

// Collect FASTQ files from sample-specific folders
FASTQ = Channel.from()
FASTQ_single = Channel.from()
fastqDirectory = file("${params.FASTQ}")
fastqDirectory.eachDir { sampleDirectory ->
	sample = sampleDirectory.name
	R1 = []
	R2 = []
	sampleDirectory.eachFileMatch(~/.*_R1.fastq.gz/) { R1_file ->
		R2_name = R1_file.name.replaceFirst(/_R1(.fastq\.gz)$/, '_R2$1')
		R2_file = file("${params.FASTQ}/${sample}/${R2_name}")
		R1.add(R1_file)
		R2.add(R2_file)
	}
	FASTQ = FASTQ.concat( Channel.from([ "R1": R1, "R2": R2, "sample": sample ]) )
	FASTQ_single = FASTQ_single.concat( Channel.from(R1) )
	FASTQ_single = FASTQ_single.concat( Channel.from(R2) )
}

// Build RG line from 1st read of each FASTQ file
process FASTQ {
	
	cpus 1
	memory '50 MB'
	scratch { scratchMode }
	
	input:
	set file(R1), file(R2), val(sample) from FASTQ
	
	output:
	set file(R1), file(R2), val(sample), stdout into FASTQ_STAR1
	set file(R1), file(R2), val(sample), stdout into FASTQ_STAR2
	
	"""
	# Get FASTQ sets from Nextflow (FIXME not space-proof)
	R1=($R1)
	R2=($R2)
	
	# Build @RG line
	RG=""
	for i in \${!R1[*]}
	do
	   # Sequencing batch ID
	   ID1=\$(gunzip -c \${R1[\$i]} | head -1 | cut -d":" -f1-4 | sed 's/^@//')
	   ID2=\$(gunzip -c \${R2[\$i]} | head -1 | cut -d":" -f1-4 | sed 's/^@//')
	   if [ \$ID1 != \$ID2 ]; then echo "Wrong R1/R2 matching" >&2; exit 1; fi
	   
	   # Mate identifier (R1/R2)
	   M1=\$(gunzip -c \${R1[\$i]} | head -1 | cut -d" " -f2 | cut -d":" -f1)
	   M2=\$(gunzip -c \${R2[\$i]} | head -1 | cut -d" " -f2 | cut -d":" -f1)
	   if [ \$M1 != "1" ]; then echo "R1 FASTQ file does not contain R1 reads" >&2; exit 1; fi
	   if [ \$M2 != "2" ]; then echo "R2 FASTQ file does not contain R2 reads" >&2; exit 1; fi
	   
	   # Barcode
	   barcode=\$(gunzip -c \${R1[\$i]} | head -1 | cut -d" " -f2 | cut -d":" -f4)
	   
	   # Complement
	   if [[ \${R1[\$i]} =~ _complement_ ]]; then suffix="_C"; else suffix=""; fi
	   
	   # Add one element to @RG, in R1 / R2 order
	   RG="\$RG , ID:${sample}\${suffix}\tBC:\${barcode}\tCN:${params.RG_CN}\tPL:${params.RG_PL}\tPM:${params.RG_PM}\tPU:\${ID1}\tSM:${sample}"
	done
	RG=\${RG:3}
	echo -n \$RG
	"""
}

// Run FastQC on individual FASTQ files
process FastQC {
	
	cpus 1
	memory '4 GB'
	time '40m'
	scratch { scratchMode }
	publishDir path: "$publishDir/QC/FastQC", mode: publishMode
	
	input:
	file FASTQ from FASTQ_single
	file adapters from file("$baseDir/in/adapters.tab")
	
	output:
	file "*.zip" into QC_FASTQC
	
	"""
	"${params.fastqc}" $FASTQ --adapters "$adapters" -o "."
	"""
}

// Build STAR index
// 2019-08-28 CALYM : 27% of 48 CPU usage, 40 GB RAM (scratch = false), 44 min
process STAR_index {
	
	cpus { params.CPU_STAR_index }
	memory '45 GB'
	time '1h'
	storeDir { storeDir }
	scratch { scratchMode }
	
	input:
	file genomeFASTA from file(params.genomeFASTA)
	file genomeGTF from file(params.genomeGTF)
	
	output:
	file("${params.genome}_raw") into rawGenome
	
	"""
	mkdir -p "./${params.genome}_raw"
	"${params.STAR}" \
		--runThreadN ${params.CPU_STAR_index} \
		--runMode genomeGenerate \
		--genomeDir "./${params.genome}_raw" \
		--genomeFastaFiles "$genomeFASTA" \
		--sjdbGTFfile "$genomeGTF" \
		--sjdbOverhang ${params.readLength - 1}
	mv "Log.out" "./${params.genome}_raw"
	"""
}

// STAR first pass
// FIXME shared-memory
process STAR_pass1 {
	
	cpus { params.CPU_STAR_align }
	memory '35 GB'
	time '1h'
	scratch { scratchMode }
	
	input:
	set file(R1), file(R2), val(sample), val(RG) from FASTQ_STAR1
	file rawGenome
	
	output:
	file("./SJ_${sample}.out.tab") into SJ
	
	"""
	mkdir -p "./$sample"
	"${params.STAR}" \
		--runThreadN ${params.CPU_STAR_align} \
		--twopassMode None \
		--genomeDir "$rawGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn "${R1.join(",")}" "${R2.join(",")}" \
		--readFilesCommand gunzip -c \
		--outFilterMultimapNmax 3 \
		--outFileNamePrefix "./" \
		--outSAMtype None
	mv ./SJ.out.tab ./SJ_${sample}.out.tab
	"""
}

// Build a new genome from STAR pass 1
process STAR_reindex {
	
	cpus 2
	memory '70 GB'
	time '1h'
	publishDir path: publishDir, mode: publishMode
	scratch { scratchMode }
	
	input:
	file SJ from SJ.collect()
	file R1 from file('in/dummy_R1.fastq')
	file R2 from file('in/dummy_R2.fastq')
	file genomeGTF from file(params.genomeGTF)
	file rawGenome
	
	output:
	file("./${params.genome}_reindexed") into reindexedGenome
	
	"""
	mkdir -p "./reindex"
	"${params.STAR}" \
	   --runThreadN 2 \
	   --genomeDir "$rawGenome" \
	   --readFilesIn "$R1" "$R2" \
	   --sjdbFileChrStartEnd $SJ \
	   --limitSjdbInsertNsj 2000000 \
	   --sjdbInsertSave All \
	   --sjdbGTFfile "$genomeGTF" \
	   --outFileNamePrefix "./reindex/" \
	   --outSAMtype None
	mv "./reindex/Log.out" "./reindex/_STARgenome/"
	mv ./reindex/_STARgenome/ ./${params.genome}_reindexed/
	"""
}

// STAR second pass
// FIXME shared-memory
process STAR_pass2 {
	
	cpus { params.CPU_STAR_align }
	memory '35 GB'
	time '3h'
	publishDir path: "$publishDir/BAM", pattern: "./${sample}.RNA.bam", mode: publishMode
	publishDir path: "$publishDir/QC/STAR", pattern: "./${sample}.log.out", mode: publishMode
	scratch { scratchMode }
	
	input:
	set file(R1), file(R2), val(sample), val(RG) from FASTQ_STAR2
	file reindexedGenome
	
	output:
	set val(sample), file("./${sample}.DNA.bam") into genomic_BAM
	set val(sample), file("./${sample}.RNA.bam") into transcriptomic_BAM
	set val(sample), file("./${sample}_SJ.out.tab") into junctions_STAR
	file "./${sample}_Log.final.out" into QC_STAR
	
	"""
	mkdir -p "./$sample"
	"${params.STAR}" \
		--runThreadN ${params.CPU_STAR_align} \
		--twopassMode None \
		--genomeDir "$reindexedGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn "${R1.join(",")}" "${R2.join(",")}" \
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
	"""
	// --chimSegmentMin ...
	// --chimOutType WithinBAM
}

// Genomically sort and index
process BAM_sort {
	
	cpus 4
	memory '4 GB'
	time '2h'
	publishDir path: "$publishDir/BAM", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM) from genomic_BAM
	
	output:
	set val(sample), file("*.sorted.bam"), file("*.sorted.bam.bai") into BAM_rnaSeqMetrics
	set val(sample), file("*.sorted.bam"), file("*.sorted.bam.bai") into BAM_featureCounts
	set val(sample), file("*.sorted.bam"), file("*.sorted.bam.bai") into BAM_secondary
//	set val(sample), file("*.sorted.bam"), file("*.sorted.bam.bai") into BAM_markDuplicates
	
	"""
	# Get file name from Nextflow
	BAM="$BAM"
	
	# Sort
	"${params.samtools}" sort -o \${BAM%.bam}.sorted.bam -T ./${sample} -@ 3 \$BAM
	
	# Index
	"${params.samtools}" index \${BAM%.bam}.sorted.bam
	"""
}

// Prepare refFlat file for Picard
// 2019-08-28 CALYM : 100% of 1 CPU usage, 378 MB RAM (scratch = false), 2.8 min
process gtfToRefFlat {
	
	cpus 1
	memory '400 MB'
	time '5m'
	storeDir { storeDir }
	scratch { scratchMode }
	
	input:
	file genomeGTF from file(params.genomeGTF)
	
	output:
	file '*.refFlat' into genomeRefFlat
	
	"""
	Rscript "${baseDir}/scripts/gtfToRefFlat.R" "$genomeGTF" "${genomeGTF}.refFlat"
	"""	
}

// Prepare rRNA interval list file for Picard
// 2019-08-28 CALYM : 96% of 1 CPU usage, 23.65 MB RAM (scratch = false), 0 min
process rRNA_interval {
	
	cpus 1
	memory '50 MB'
	time '1m'
	storeDir { storeDir }
	scratch { scratchMode }
	
	input:
	file genomeGTF from file(params.genomeGTF)
	file rawGenome
	
	output:
	file './rRNA.interval_list' into rRNA_interval
	
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
	memory '12 GB'
	time '90m'
	publishDir path: "$publishDir/QC/rnaSeqMetrics", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM), file(BAI) from BAM_rnaSeqMetrics
	file rRNA_interval
	file genomeRefFlat
	
	output:
	file "./${sample}.RNA_Metrics" into QC_rnaSeqMetrics
	
	"""
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "${params.picard}" CollectRnaSeqMetrics \
		INPUT=$BAM \
		OUTPUT="./${sample}.RNA_Metrics" \
		REF_FLAT="$genomeRefFlat" \
		RIBOSOMAL_INTERVALS="$rRNA_interval" \
		STRAND_SPECIFICITY="${params.strandness}" \
		ASSUME_SORTED=true
	"""
}

// Count reads in transcripts using featureCounts
process featureCounts {
	
	cpus 2
	memory '800 MB'
	time '1h'
	publishDir path: "$publishDir/featureCounts", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM), file(BAI) from BAM_featureCounts
	file genomeGTF from file(params.genomeGTF)
	file genomeRefFlat
	
	output:
	file "./annotation.csv" into featureCounts_annotation
	file "./${sample}_counts.rds" into featureCounts_counts
	file "./${sample}_stats.rds" into featureCounts_stats
	
	"""
	#!/usr/bin/env Rscript
	
	# Dependency
	.libPaths("${params.R_packages}")
	library(Rsubread)
	
	# Count reads in genes (~5 minutes with 8 threads)
	dir.create("./tmp")
	out <- featureCounts(
		files = "$BAM",
		annot.ext = "$genomeGTF",
		isGTFAnnotationFile = TRUE,
		allowMultiOverlap = FALSE,
		minMQS = 0L,
		strandSpecific = ${params.strandSpecificity},
		isPairedEnd = TRUE,
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
	memory '800 MB'
	time '10m'
	publishDir path: "$publishDir", pattern: "all_counts.rds", mode: publishMode
	publishDir path: "$publishDir/QC", pattern: "*.yaml", mode: publishMode
	scratch { scratchMode }
	
	input:
	file 'annotation' from featureCounts_annotation.first()
	file 'countFiles' from featureCounts_counts.collect()
	
	output:
	file 'all_counts.rds' into countMatrix
	file './edgeR.yaml' into QC_edgeR_general
	file './edgeR_mqc.yaml' into QC_edgeR_section
	
	"""
	Rscript "${baseDir}/scripts/edgeR.R" "$annotation" "." "${params.R_packages}" $countFiles
	"""	
}

// Estimate insert size distribution
process insertSize {
	
	cpus 2
	memory '500 MB'
	time '10m'
	publishDir path: "$publishDir/QC/insertSize", pattern: "./${sample}_mqc.yaml", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM) from transcriptomic_BAM
	
	output:
	file "./${sample}_mqc.yaml" into QC_insert
	
	"""
	Rscript "${baseDir}/scripts/insertSize.R" "$sample" "$BAM" "${params.samtools}" > "./${sample}_mqc.yaml"
	"""	
}

/*
// Picard MarkDuplicates (remove optical duplicates but keep and mark library-related duplicates)
// FIXME use as many CPUs as available, whatever the options -> disabled, not used so far anyway
// FIXME sort and index after markDuplicates, not before
process markDuplicates {
	
	cpus 1
	publishDir path: "$publishDir/QC/markDuplicates", pattern: "./${sample}.txt", mode: publishMode
	publishDir path: "$publishDir/BAM", pattern: "./${BAM.getBaseName()}.MD.bam", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM), file(BAI) from BAM_markDuplicates
	
	output:
	file "./${sample}.txt" into QC_markDuplicates
	file "./${BAM.getBaseName()}.MD.bam" into BAM_marked
	
	"""
	java -Xmx4G -Duser.country=US -Duser.language=en -jar "${params.picard}" MarkDuplicates \
		TMP_DIR="." \
		INPUT=$BAM \
		OUTPUT="./${BAM.getBaseName()}.MD.bam" \
		METRICS_FILE="./${sample}.txt" \
		ASSUME_SORT_ORDER="coordinate" \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_SEQUENCING_DUPLICATES="true" \
		REMOVE_DUPLICATES="false" \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=50
	"""
}
*/

// Quantify secondary alignments with SAMtools
// TODO : general stats
process secondary {
	
	cpus 1
	memory '200 MB'
	time '1h'
	publishDir path: "$publishDir/QC/secondary", pattern: "./${sample}_mqc.yaml", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(BAM), file(BAI) from BAM_secondary
	
	output:
	file "./${sample}_mqc.yaml" into QC_secondary
	
	"""
	# Total read count (fast, from index)
	total=\$(${params.samtools} idxstats $BAM | awk 'BEGIN{x=0} {x=x+\$3+\$4} END{print x}')
	
	# Secondary alignments
	sec=\$(${params.samtools} view -c -f 0x100 $BAM)
	
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
	memory '4 GB'
	time '10m'
	publishDir path: "$publishDir/QC", mode: publishMode
	scratch { scratchMode }
	
	input:
	file conf from file('in/multiqc.conf')
	file 'edgeR.yaml' from QC_edgeR_general
	file 'edgeR_mqc.yaml' from QC_edgeR_section
	file 'STAR/*' from QC_STAR.collect()
	file 'FASTQC/*' from QC_FASTQC.collect()
//	file 'markDuplicates/*' from QC_markDuplicates.collect()
	file 'rnaSeqMetrics/*' from QC_rnaSeqMetrics.collect()
	file 'insertSize/*' from QC_insert.collect()
	file 'secondary/*' from QC_secondary.collect()
	
	output:
	file "*_multiqc_report_data.zip" into MultiQC_data
	file "*_multiqc_report.html" into MultiQC_report
	
	"""
	multiqc --title "${params.MQC_title}" --comment "${params.MQC_comment}" --outdir "." --config "${conf}" --config "./edgeR.yaml" --zip-data-dir --interactive --force "."
	"""
}

// Reshape STAR junction file into a Rgb table
process junctions {
	
	cpus 1
	memory '300 MB'
	time '10m'
	publishDir path: "$publishDir/junctions", pattern: "./${sample}.rdt", mode: publishMode
	scratch { scratchMode }
	
	input:
	set val(sample), file(SJ_tab) from junctions_STAR
	
	output:
	file "./${sample}.rdt" into junctions_Rgb
	
	"""
	#!/usr/bin/env Rscript
	
	# Dependency
	.libPaths("${params.R_packages}")
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
