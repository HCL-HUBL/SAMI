process splicing_filter {
	cpus 10
	time { 30.minute * task.attempt }
	memory { 10.GB + 5.GB * task.attempt }
	publishDir "${params.out}/splicing", mode: params.publish

	input:
	path(exons)
	path('*')
	path('*')
	path('*')
	path(targetGTF)
	val(dir)

	output:
	path("${dir}"), emit: splicing_output

	"""
	Rscript --vanilla "${projectDir}/scripts/splicing_filter.R" ${cpus} "$targetGTF" "$exons" ${params.plot} ${params.fusions} ${params.min_I} ${params.min_PSI} "$params.symbols" "$params.classes" "$params.focus"
	"""
}
