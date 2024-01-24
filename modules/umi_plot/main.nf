process umi_plot {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    tuple(val(sample), file(umiHist))

    output:
    path(env(outQC)), emit: QC_umi

    """
    Rscript --vanilla "${baseDir}/scripts/umi_stat.R" "$sample" "$umiHist"

    outQC="${sample}_mqc.yaml"
    """
}
