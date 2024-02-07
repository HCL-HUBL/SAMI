process duplication_umi_based {
    cpus 1
    label 'monocore'
    label 'retriable'

    input:
    path('*')
    path('*')

    output:
    path("duplication_umi.yaml"), emit: dup_umi

    """
    bash "${projectDir}/scripts/duplication_umi.sh" "${params.out}"
    """
}
