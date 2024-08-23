process duplication_umi_based {
	cpus 1
	time { 30.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	path('*')
	path('*')

	output:
	path("duplication_umi.yaml"), emit: YAML

	shell:
	template 'duplication_umi.sh'
}
