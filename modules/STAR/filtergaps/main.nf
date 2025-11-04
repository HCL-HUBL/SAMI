process star_filtergaps {
	tag "$sample"
	
	cpus 1
	time { 30.minute * task.attempt }
	memory { 2.GB * task.attempt }

	input:
	tuple val(sample), path(junctions), path(depth), path(index)
	val(threshold)
	
	output:
	path("out/*"), emit: junctions

	shell:
	template 'filtergaps.R'
}
