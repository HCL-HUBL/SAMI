process edgeR {
	cpus 1
	time { 15.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	path(annotation)
	path(countFiles)

	output:
	tuple path("counts.tsv"), path("CPM.tsv"), path("RPK.tsv"), path(annotation), emit: TSV
	path("edgeR.yaml"), emit: YAML_general
	path("edgeR_mqc.yaml"), emit: YAML_section

	shell:
	template 'edgeR.R'
}
