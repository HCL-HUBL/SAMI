process annotation {

    cpus 1
    label 'monocore'
    label 'retriable'

    when:
    params.splicing

    input:
    path(genomeGTF)

    output:
    path("${genomeGTF}.introns.rds"), emit: introns
    path("${genomeGTF}.exons.rdt"), emit: exons_collect
    path("${genomeGTF}.exons.rdt"), emit: exons_filter
    path("${genomeGTF}.genes.rdt"), emit: genes

    """
    Rscript --vanilla "${projectDir}/scripts/annotation.R" "$genomeGTF" "$params.species" "$params.genome" "$params.chromosomes"
    """
}
