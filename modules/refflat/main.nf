process refflat {

    cpus 1
    label 'monocore'
    label 'nonRetriable'

    input:
    path(GTF)

    output:
    path("${GTF}.refFlat"), emit: refFlats

    """
    Rscript --vanilla "${projectDir}/scripts/gtfToRefFlat.R" "$GTF" "${GTF}.refFlat"
    """
}
