process refflat {
    cpus 1
    label 'monocore'
    label 'nonRetriable'
    storeDir params.store

    input:
    path(GTF)

    output:
    path("${GTF}.refFlat"), emit: refFlats

    """
    Rscript --vanilla "${projectDir}/scripts/gtfToRefFlat.R" "$GTF" "${GTF}.refFlat"
    """
}
