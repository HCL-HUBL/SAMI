process annotation {
    cpus 1
    label 'monocore'
    label 'retriable'
    storeDir params.store

    when:
    params.splicing

    input:
    path(genomeGTF)

    output:
    path("${genomeGTF}.introns.rds"), emit: introns
    path("${genomeGTF}.exons.rdt"), emit: exons
    path("${genomeGTF}.genes.rdt"), emit: genes

    """
    Rscript --vanilla "${projectDir}/scripts/annotation.R" "$genomeGTF" "$params.species" "$params.genome" "$params.chromosomes"
    """
}
