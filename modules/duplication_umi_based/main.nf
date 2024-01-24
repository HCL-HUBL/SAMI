process duplication_umi_based {

    cpus 1
    label 'monocore'
    label 'retriable'
    publishDir "${params.out}/QC", mode: "copy"

    input:
    path('*')
    path('*')

    output:
    path(env(outYAML)), emit: dup_umi

    """
    bash "${baseDir}/scripts/duplication_umi.sh" "${params.out}"
    outYAML="duplication_umi.yaml"
    """
}
