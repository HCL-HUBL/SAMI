process featurecounts {
	tag "$sample"

	cpus 2
	time { 1.hour * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.out}/expression", mode: params.publish, pattern: "annotation.tsv"

	input:
	tuple val(sample), val(type), path(BAM), path(BAI)
	path(targetGTF)

	output:
	path("annotation.tsv"), emit: annotation
	path("${sample}_counts.rds"), emit: counts

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
