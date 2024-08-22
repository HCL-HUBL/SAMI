process markduplicates {
    tag "$sample"

    cpus 1
	time { 15.minute * task.attempt }
	memory { 5.GB * task.attempt }

    input:
    tuple val(sample), val(type), path(BAM)

    output:
    path("${sample}.txt"), emit: txt
    tuple val(sample), val(type), path("${BAM.getBaseName()}.MD.bam"), emit: BAM

    """
    java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" MarkDuplicates \
        --TMP_DIR "." \
        --INPUT "$BAM" \
        --OUTPUT "${BAM.getBaseName()}.MD.bam" \
        --METRICS_FILE "${sample}.txt" \
        --ASSUME_SORT_ORDER "queryname" \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
        --REMOVE_SEQUENCING_DUPLICATES "false" \
        --REMOVE_DUPLICATES "false" \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 50 \
        --PROGRAM_RECORD_ID null
    """
}
