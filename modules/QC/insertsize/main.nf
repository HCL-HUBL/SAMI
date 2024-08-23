process insertsize {
	tag "$sample"

	cpus 1
	time { 10.minute * task.attempt }
	memory { 500.MB * task.attempt }

	when:
	type == "paired"

	input:
	tuple val(sample), val(type), path(isize)

	output:
	path("${sample}_mqc.yaml"), emit: YAML

	shell:
	template 'insertSize.R'
}
