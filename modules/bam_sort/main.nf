process bam_sort {
    tag "$sample"

    cpus 4
    label 'multicore'
    label 'retriable'
    publishDir "${params.out}/bam_splicing", mode: ${params.publish}

    input:
    tuple val(sample), val(type), path(BAM)

    output:
    tuple val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai"), emit: BAM_sorted
    path("${BAM.getBaseName()}.sort.bam"), emit: onlyBAM_sorted
    path("${BAM.getBaseName()}.sort.bai"), emit: BAI_splicing

    """
    # Sort
    samtools sort -o "${BAM.getBaseName()}.sort.bam" -T ./${sample} -@ 3 "$BAM"

    # Index
    samtools index "${BAM.getBaseName()}.sort.bam"
    mv "${BAM.getBaseName()}.sort.bam.bai" "${BAM.getBaseName()}.sort.bai"
    """
}
