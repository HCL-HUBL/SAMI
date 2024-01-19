process umi_plot {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    tuple(val(sample), file(umiHist)) from UMI_stat

    output:
    path(env(outQC)), emit: QC_umi

    if(params.umi) {
        """
        Rscript --vanilla "${baseDir}/scripts/umi_stat.R" "$sample" "$umiHist"

        outQC="${sample}_mqc.yaml"
        """
    } else {
        """
        outQC="${baseDir}/in/dummy.tsv"
        """
    }
}
