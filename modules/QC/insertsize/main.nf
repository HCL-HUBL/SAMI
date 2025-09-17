process insertsize {
	tag "$sample"

	cpus 1
	time { 15.minute * task.attempt }
	memory { 2.GB * task.attempt }

	when:
	type == "paired"

	input:
	tuple val(sample), val(type), path(isize)
	path("*")

	output:
	path("${sample}_mqc.yaml"), emit: YAML
	path("isize_table.yaml"), emit: TABLE

	shell:
	template 'insertSize.R'
}
