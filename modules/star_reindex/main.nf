process star_reindex {

    cpus 2
    label 'multicore'
    label 'retriable'
    publishDir params.out, mode: "copy"

    input:
    path(SJ)
    path(genomeGTF)
    path(rawGenome)

    output:
    path("${params.genome}_${params.title}"), emit: reindexedGenome

    """
    mkdir -p "./reindex"
    STAR \
        --runThreadN 2 \
        --genomeDir "$rawGenome" \
        --readFilesIn "${baseDir}/in/dummy_R1.fastq" "${baseDir}/in/dummy_R2.fastq" \
        --sjdbFileChrStartEnd $SJ \
        --limitSjdbInsertNsj 5000000 \
        --sjdbInsertSave All \
        --sjdbGTFfile "$genomeGTF" \
        --outFileNamePrefix "./reindex/" \
        --outSAMtype None
    mv "./reindex/Log.out" "./reindex/_STARgenome/"
    mv "./reindex/_STARgenome/" "./${params.genome}_${params.title}/"
    """
}
