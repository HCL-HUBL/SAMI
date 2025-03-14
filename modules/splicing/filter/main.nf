process splicing_filter {
	cpus 10
	time { 30.minute * task.attempt }
	memory { 10.GB + 5.GB * task.attempt }
	publishDir "${params.output}/splicing", mode: params.publish

	input:
	path(exons)
	tuple path("I.rds"), path("S.rds"), path("groups.rds"), path("sites.rds"), path("events.rds")
	path("depth.bed")
	path(targetGTF)
	val(dir)
	val(plot)
	val(fusions)
	val(min_I)
	val(min_PSI)
	val(symbols)
	val(classes)

	output:
	path("${dir}"), emit: dir

	shell:
	template 'splicing_filter.R'
}
