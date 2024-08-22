process umi_table {

    cpus 1
    label 'retriable'

    input:
    path("*")

    output:
    path("umi_table_mqc.yaml"), emit: YAML

    """
    Rscript --vanilla "${projectDir}/scripts/umi_table.R"
    """
}
