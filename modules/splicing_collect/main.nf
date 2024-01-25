process splicing_collect {

    cpus { params.CPU_splicing }
    label 'multicore'
    label 'nonRetriable'
    publishDir "${params.out}/splicing", mode: "copy"

    when:
    params.splicing

    input:
    path(genes)
    path(exons)
    path(introns)
    path('junctionFiles/*')
    path('chimericFiles/*')
    path("transcripts.tsv")

    output:
    path("*.rds"), emit: splicing_events

    """
    Rscript --vanilla "${projectDir}/scripts/splicing_collect.R" ${params.CPU_splicing} "$genes" "$exons" "$introns" "$params.chromosomes" $params.min_reads_unknown "transcripts.tsv"
    """
}
