process insertsize_table {
	cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	path("*")

	output:
	path("isize_table.yaml"), emit: YAML

	shell:
	template 'insert_table.R'
}
