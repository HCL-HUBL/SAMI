process umi_plot {
    tag "$sample"

    cpus 1
    label 'retriable'

    input:
    tuple val(sample), path(umiHist)

    output:
    path("${sample}_mqc.yaml"), emit: YAML

    """
    Rscript --vanilla "${projectDir}/scripts/umi_stat.R" "$sample" "$umiHist"
    """
}
