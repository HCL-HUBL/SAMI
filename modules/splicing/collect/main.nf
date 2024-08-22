process splicing_collect {
	cpus 10
	time { 10.minute * task.attempt }
	memory { 10.GB + 5.GB * task.attempt }

	input:
	path(genes)
	path(exons)
	path(introns)
	path('junctionFiles/*')
	path('chimericFiles/*')
	path("transcripts.tsv")
	val(chromosomes)
	val(min_reads_unknown)
	
	output:
	path("*.rds"), emit: RDS

	"""
	Rscript --vanilla "${projectDir}/scripts/splicing_collect.R" ${task.cpus} "$genes" "$exons" "$introns" "$chromosomes" $min_reads_unknown "transcripts.tsv"
	"""
}
