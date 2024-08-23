process versions {
	cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	val gitVersion

	output:
	path("SAMI_mqc_versions.yaml"), emit: YAML
	
	shell:
	template 'versions.bash'
}
