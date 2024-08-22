process multiqc {
    cpus 1
	time { 20.minute * task.attempt }
	memory { 4.GB * task.attempt }
    publishDir "${params.out}/QC", mode: params.publish

    input:
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
    path("${params.MQC_title}_multiqc_report_data.zip")
    path("${params.MQC_title}_multiqc_report.html")

    """
    multiqc --title "${params.MQC_title}" --comment "${params.MQC_comment}" --outdir "." --config "${projectDir}/in/multiqc.conf" --config "./edgeR.yaml" --config "./umi_table_mqc.yaml" --config "./duplication_umi.yaml" --config "./isize_table_mqc.yaml" --zip-data-dir --interactive --force "."
    """
}
