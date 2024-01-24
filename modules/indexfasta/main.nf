process indexfasta {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir params.store, mode: "copy"
    scratch params.scratch

    input:
    path(genomeFASTA)

    output:
    tuple path(genomeFASTA), path("${genomeFASTA.getBaseName()}.dict"), path("${genomeFASTA}.fai"), emit: indexedFASTA_splitN
    tuple path(genomeFASTA), path("${genomeFASTA.getBaseName()}.dict"), path("${genomeFASTA}.fai"), emit: indexedFASTA_BQSR
    tuple path(genomeFASTA), path("${genomeFASTA.getBaseName()}.dict"), path("${genomeFASTA}.fai"), emit: indexedFASTA_Mutect2
    tuple path(genomeFASTA), path("${genomeFASTA.getBaseName()}.dict"), path("${genomeFASTA}.fai"), emit: indexedFASTA_MergeBamAlignment

    """
    # Dictionnary
    java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CreateSequenceDictionary \
        REFERENCE="$genomeFASTA" \
        OUTPUT="${genomeFASTA.getBaseName()}.dict"
    # Index
    samtools faidx "$genomeFASTA"
    """
}
