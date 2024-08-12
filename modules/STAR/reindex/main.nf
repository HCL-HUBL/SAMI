process star_reindex {

    cpus 2
    label 'multicore'
    label 'retriable'

    input:
    path(SJ)
    path(rawGenome)
    path(genomeGTF)

    output:
    path("${params.genome}_${params.title}"), emit: genome

    """
    mkdir -p "./reindex"
    STAR \
        --runThreadN 2 \
        --genomeDir "$rawGenome" \
        --readFilesIn "${projectDir}/in/dummy_R1.fastq" "${projectDir}/in/dummy_R2.fastq" \
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
