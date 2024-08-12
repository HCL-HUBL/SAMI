process star_index {
    cpus { params.CPU_index }
    label 'multicore'
    label 'nonRetriable'
    storeDir params.store

    input:
    path(genomeFASTA)
    path(genomeGTF)

    output:
    path("${params.genome}_raw"), emit: genome
    path("${params.genome}_raw/chrNameLength.txt"), emit: chrom

    """
    mkdir -p "./${params.genome}_raw"
    STAR \
        --runThreadN ${params.CPU_index} \
        --runMode genomeGenerate \
        --genomeDir "./${params.genome}_raw" \
        --genomeFastaFiles "$genomeFASTA" \
        --sjdbGTFfile "$genomeGTF"
    """
}
