process softclipping {
	tag "$sample"

	cpus 1
	time { 30.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	tuple val(sample), val(type), path(BAM), path(BAI)

	output:
	path("${sample}_*_mqc.yaml"), emit: YAML

	shell:
	template 'softClipping.R'
}
