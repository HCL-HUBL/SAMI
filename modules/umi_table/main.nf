process umi_table {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    path("*")

    output:
    path(env(outYAML)), emit: QC_umi_table

    if(params.umi) {
        """
        Rscript --vanilla "${baseDir}/scripts/umi_table.R"
        outYAML=umi_table_mqc.yaml
        """
    } else {
        """
        outYAML="${baseDir}/in/dummy.tsv"
        """
    }
}
