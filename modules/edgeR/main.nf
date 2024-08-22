process edgeR {
    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/expression", mode: params.publish, pattern: "*.tsv"

    input:
    path(annotation)
    path(countFiles)

    output:
    path("counts.tsv"), emit: counts
    path("CPM.tsv"), emit: CPM
    path("RPK.tsv"), emit: RPK
    path("edgeR.yaml"), emit: YAML_general
    path("edgeR_mqc.yaml"), emit: YAML_section

    """
    Rscript --vanilla "${projectDir}/scripts/edgeR.R" "$annotation" "." $countFiles
    """
}
