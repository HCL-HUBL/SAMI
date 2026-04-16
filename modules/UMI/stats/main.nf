process umi_stats {
	cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	path('histograms/*')

	output:
	path("umi-table.yaml"), emit: table
	path("umi-plot_mqc.yaml"), emit: plot

	shell:
	template 'umi_stats.R'
}
