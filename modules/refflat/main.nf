process refflat {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir params.store, mode: "copy"

    input:
    path(GTF)

    output:
    path("${GTF}.refFlat"), emit: refFlats

    """
    Rscript --vanilla "${projectDir}/scripts/gtfToRefFlat.R" "$GTF" "${GTF}.refFlat"
    """
}
