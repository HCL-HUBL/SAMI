process bam_sort {

    cpus 4
    label 'multicore'
    label 'retriable'
    storeDir { "${params.out}/BAM" }

    input:
    tuple(val(sample), val(type), path(BAM))

    output:
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai")), emit: BAM_sorted
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai")), emit: BAM_rnaSeqMetrics
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai")), emit: BAM_featureCounts
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai")), emit: BAM_secondary
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai")), emit: BAM_softClipping
    path("${BAM.getBaseName()}.sort.bam"), emit: BAM_splicing
    path("${BAM.getBaseName()}.sort.bam"), emit: BAM_dup2
    path("${BAM.getBaseName()}.sort.bai"), emit: BAI_splicing
    path("${BAM.getBaseName()}.sort.clean"), emit: BAM_sort_clean

    """
    # Sort
    samtools sort -o "${BAM.getBaseName()}.sort.bam" -T ./${sample} -@ 3 "$BAM"

    # Index
    samtools index "${BAM.getBaseName()}.sort.bam"
    mv "${BAM.getBaseName()}.sort.bam.bai" "${BAM.getBaseName()}.sort.bai"

    # Link input BAM for cleaning
    ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.sort.clean"
    """
}
