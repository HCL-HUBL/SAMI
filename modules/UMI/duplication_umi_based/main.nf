process duplication_umi_based {
    cpus 1
    time { 30.minute * task.attempt }
	memory { 500.MB * task.attempt }

    input:
    path('*')
    path('*')

    output:
    path("duplication_umi.yaml"), emit: YAML

    """
    bash "${projectDir}/scripts/duplication_umi.sh" "${params.out}"
    """
}
