process umi_plot {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    tuple val(sample), path(umiHist)

    output:
    path(outQC), emit: QC_umi

    """
    Rscript --vanilla "${projectDir}/scripts/umi_stat.R" "$sample" "$umiHist"

    outQC="${sample}_mqc.yaml"
    """
}
