process insertsize {

    cpus 1
    label 'monocore'
    label 'retriable'
    publishDir "${params.out}/QC/insertSize", mode: "copy"

    when:
    type == "paired"

    input:
    tuple val(sample), val(type), path(isize)

    output:
    path("${sample}_mqc.yaml"), emit: QC_insert

    """
    Rscript --vanilla "${baseDir}/scripts/insertSize.R" "$sample" "$isize" > "./${sample}_mqc.yaml"
    """
}
