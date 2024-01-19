process versions {

	cpus 1
	label 'nonRetriable'
	label 'monocore'
	executor 'local'

	publishDir "${params.out}/QC", mode: "copy"

	input:
	path(script) from file("${baseDir}/scripts/versions.bash")

	output:
	path("SAMI_mqc_versions.yaml"), emit: versions

	"""
	#!/bin/bash
	rm -f "SAMI_mqc_versions.yaml"
	echo "Command: '${workflow.commandLine}'" >> "SAMI_mqc_versions.yaml"
	echo "Container: '${workflow.container}'" >> "SAMI_mqc_versions.yaml"
	echo "Nextflow: '${nextflow.version}'" >> "SAMI_mqc_versions.yaml"
	echo "SAMI: '$gitVersion'" >> "SAMI_mqc_versions.yaml"
	bash $script >> "SAMI_mqc_versions.yaml"
	"""
}
