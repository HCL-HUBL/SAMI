process splicing_harvest {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	tuple val(sample), val(type), path(BAM), path(BAI)
	tuple path(genomeFASTA), path(genomeDICT), path(genomeFAI)
	val(flag_exclude)
	
	output:
	path("${BAM}.tsv"), emit: TSV

	"""
	harvest "${BAM}" "${genomeFASTA}" 0 "${flag_exclude}" > "${BAM}.tsv"
	"""
}
