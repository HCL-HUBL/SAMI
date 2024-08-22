process umi_plot {
    tag "$sample"

    cpus 1
    time { 5.minute * task.attempt }
	memory { 500.MB * task.attempt }

    input:
    tuple val(sample), path(umiHist)

    output:
    path("${sample}_mqc.yaml"), emit: YAML

    """
    Rscript --vanilla "${projectDir}/scripts/umi_stat.R" "$sample" "$umiHist"
    """
}
