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
	tuple path("I.rds"), path("S.rds"), path("groups.rds"), path("sites.rds"), path("events.rds"), path("splicing.rds"), emit: RDS
	path("depth.bed"), emit: BED

	shell:
	template 'splicing_aggregate.R'
}
