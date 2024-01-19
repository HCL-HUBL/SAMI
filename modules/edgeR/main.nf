process edgeR {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/edgeR", mode: "copy"

    when:
    params.finalize

    input:
    path(annotation)
    path(countFiles)

    output:
    path("counts.tsv"), emit: counts
    path("CPM.tsv"), emit: CPM
    path("RPK.tsv"), emit: RPK
    path("edgeR.yaml"), emit: QC_edgeR_general
    path("edgeR_mqc.yaml"), emit: QC_edgeR_section

    """
    Rscript --vanilla "${baseDir}/scripts/edgeR.R" "$annotation" "." $countFiles
    """
}
