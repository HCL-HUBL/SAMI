process umi_plot {
    tag "$sample"

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    tuple val(sample), path(umiHist)

    output:
    path("${sample}_mqc.yaml"), emit: QC_umi

    """
    Rscript --vanilla "${projectDir}/scripts/umi_stat.R" "$sample" "$umiHist"
    """
}
