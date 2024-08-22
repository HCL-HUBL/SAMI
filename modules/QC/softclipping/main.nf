process softclipping {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'retriable'

    input:
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    path("${sample}_*_mqc.yaml"), emit: YAML

    """
    Rscript --vanilla "${projectDir}/scripts/softClipping.R" "$sample" "$BAM"
    """
}
