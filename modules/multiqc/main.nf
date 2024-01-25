process multiqc {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/QC", mode: "copy"

    when:
    params.finalize

    input:
    path('edgeR.yaml')
    path('edgeR_mqc.yaml')
    path('STAR_pass1/*')
    path('STAR_pass2/*')
    path('FastQC/raw/*_fastqc.zip')
    path('FastQC/trimmed/*_fastqc.zip')
    path('markDuplicates/*')
    path('rnaSeqMetrics/*')
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
    path("${params.MQC_title}_multiqc_report_data.zip"), emit: MultiQC_data
    path("${params.MQC_title}_multiqc_report.html"), emit: MultiQC_report

    """
    multiqc --title "${params.MQC_title}" --comment "${params.MQC_comment}" --outdir "." --config "${projectDir}/in/multiqc.conf" --config "./edgeR.yaml" --config "./umi_table_mqc.yaml" --config "./duplication_umi.yaml" --config "./isize_table_mqc.yaml" --zip-data-dir --interactive --force "."
    """
}
