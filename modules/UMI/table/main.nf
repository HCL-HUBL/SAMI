process umi_table {
	cpus 1
	time { 5.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	path("*")

	output:
	path("umi_table_mqc.yaml"), emit: YAML

	shell:
	template 'umi_table.R'
}
