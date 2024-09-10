process splicing_aggregate {
	cpus 10
	time { 10.minute * task.attempt }
	memory { 10.GB + 5.GB * task.attempt }

	input:
	path(genes)
	path(exons)
	path(introns)
	path('gapFiles/*')
	path('chimericFiles/*')
	path("transcripts.tsv")
	val(chromosomes)
	val(min_reads_unknown)
	val(stranded)
	
	output:
	path("*.rds"), emit: RDS

	shell:
	template 'splicing_aggregate.R'
}
