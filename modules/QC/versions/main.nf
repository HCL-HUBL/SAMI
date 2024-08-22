process versions {
	cpus 1
	time { 1.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	val gitVersion

	output:
	path("SAMI_mqc_versions.yaml"), emit: YAML

	"""
	#!/bin/bash
	rm -f "SAMI_mqc_versions.yaml"
	echo "Command: ${workflow.commandLine}" >> "SAMI_mqc_versions.yaml"
	echo "Container: '${workflow.container}'" >> "SAMI_mqc_versions.yaml"
	echo "Nextflow: '${nextflow.version}'" >> "SAMI_mqc_versions.yaml"
	echo "SAMI: '${gitVersion}'" >> "SAMI_mqc_versions.yaml"
	bash "${projectDir}/scripts/versions.bash" >> "SAMI_mqc_versions.yaml"
	"""
}
