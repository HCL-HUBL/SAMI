process depth {
	tag "$sample"

	cpus 1
	time { 1.hour + 1.hour * task.attempt }
	memory { 1.GB * task.attempt }
	
	input:
	tuple val(sample), val(type), path(BAM), path(BAI)

	output:
	tuple val(sample), path("${BAM.getBaseName()}.tsv.gz"), path("${BAM.getBaseName()}.tsv.gz.tbi"), emit: TSV

	"""
	# Depth along whole genome
	samtools depth "$BAM" -o "${BAM.getBaseName()}.tsv"
	
	# Compress
	bgzip "${BAM.getBaseName()}.tsv"
	
	# Index
	tabix -s 1 -b 2 -e 2 "${BAM.getBaseName()}.tsv.gz"
	"""
}
