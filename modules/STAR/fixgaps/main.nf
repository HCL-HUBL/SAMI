process star_fixgaps {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	path(exons)
	tuple path(fasta), path(dict), path(fai)
	path('gapFiles/*')
	val(range)
	
	output:
	path("fixed.out.tab"), emit: junctions

	shell:
	template 'fixgaps.R'
}
