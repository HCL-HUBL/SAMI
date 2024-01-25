process umi_table {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/umi", mode: "copy"

    input:
    path("*")

    output:
    path("umi_table_mqc.yaml"), emit: QC_umi_table

    """
    Rscript --vanilla "${projectDir}/scripts/umi_table.R"
    """
}
