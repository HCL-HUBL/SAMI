process star_index {
    cpus { params.CPU_index }
    label 'multicore'
    label 'nonRetriable'
    publishDir params.store, mode: "copy"

    input:
    path(genomeFASTA)
    path(genomeGTF)

    output:
    path("${params.genome}_raw"), emit: rawGenome_pass1
    path("${params.genome}_raw"), emit: rawGenome_reindex
    path("${params.genome}_raw/chrNameLength.txt"), emit: rawGenome_chrom

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
