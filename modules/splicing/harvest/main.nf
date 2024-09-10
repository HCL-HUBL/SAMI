process splicing_harvest {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	tuple val(sample), val(type), path(BAM), path(BAI)
	
	output:
	path("${BAM}.tsv"), emit: TSV

	"""
	harvest "${BAM}" > "${BAM}.tsv"
	"""
}
