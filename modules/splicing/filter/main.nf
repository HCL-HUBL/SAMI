process splicing_filter {
	cpus 10
	time { 30.minute * task.attempt }
	memory { 10.GB + 5.GB * task.attempt }
	publishDir "${params.out}/splicing", mode: params.publish

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
	val(focus)

	output:
	path("${dir}"), emit: splicing_output

	shell:
	template 'splicing_filter.R'
}
