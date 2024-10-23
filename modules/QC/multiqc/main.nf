process multiqc {
	cpus 1
	time { 20.minute * task.attempt }
	memory { 4.GB * task.attempt }
	publishDir "${params.output}/QC", mode: params.publish

	input:
	val(title)
	val(comment)
	path('multiqc.conf')
	path('edgeR.yaml')
	path('edgeR_mqc.yaml')
	path('STAR_pass1/*')
	path('STAR_pass2/*')
	path('FastQC/raw/*_fastqc.zip')
	path('FastQC/trimmed/*_fastqc.zip')
	path('markDuplicates/*')
	path('rnaSeqMetrics/genome/*')
	path('rnaSeqMetrics/target/*')
	path('insertSize/*')
	path('secondary/*')
	path('softClipping/*')
	path('umi/*_mqc.yaml')
	path('umi_table_mqc.yaml')
	path('isize_table_mqc.yaml')
	path('cutadapt/*')
	path('duplication_umi.yaml')
	path('SAMI_mqc_versions.yaml')

	output:
	path("${title}_multiqc_report_data.zip"), emit: ZIP
	path("${title}_multiqc_report.html"), emit: HTML

	"""
	# Mandatory config files
	args=""
	args="\$args --config multiqc.conf"
	args="\$args --config edgeR.yaml"
	args="\$args --config isize_table_mqc.yaml"
	
	# Optional config files
	if [ -f "umi_table_mqc.yaml" ];   then args="\$args --config umi_table_mqc.yaml";   fi
	if [ -f "duplication_umi.yaml" ]; then args="\$args --config duplication_umi.yaml"; fi
	
	# Main arguments
	args="\$args --zip-data-dir --interactive --force ."
	
	# Call
	multiqc --title "${title}" --comment "${comment}" --outdir . \$args
	"""
}
