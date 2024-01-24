process softclipping {

    cpus 1
    label 'monocore'
    label 'retriable'
    publishDir "${params.out}/QC/softClipping", mode: "copy"

    input:
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    path("${sample}_*_mqc.yaml"), emit: QC_softClipping

    """
    Rscript --vanilla "${baseDir}/scripts/softClipping.R" "$sample" "$BAM"
    """
}
