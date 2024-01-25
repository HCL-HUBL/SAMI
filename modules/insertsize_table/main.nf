process insertsize_table {

    cpus 1
    label 'retriable'
    publishDir "${params.out}/QC/insertSize", mode: "copy"

    input:
    path("*")

    output:
    path("isize_table_mqc.yaml"), emit: median_isize_table

    """
    Rscript --vanilla "${projectDir}/scripts/insert_table.R"
    """
}
