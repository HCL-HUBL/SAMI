process insertsize {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'retriable'

    when:
    type == "paired"

    input:
    tuple val(sample), val(type), path(isize)

    output:
    path("${sample}_mqc.yaml"), emit: YAML

    """
    Rscript --vanilla "${projectDir}/scripts/insertSize.R" "$sample" "$isize" > "./${sample}_mqc.yaml"
    """
}
