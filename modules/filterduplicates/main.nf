process filterduplicates {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/filterDuplicates", mode: "copy"
    scratch { params.scratch }

    when:
    params.varcall

    input:
    tuple(val(sample), val(type), path(BAM), path(BAI))

    output:
    tuple(val(sample), val(type), path("${BAM.getBaseName()}.filter.bam"), path("${BAM.getBaseName()}.filter.bai")), emit: BAM_filtered

    """
    # Filter
    samtools view -b -F 0x400 "$BAM" -o "${BAM.getBaseName()}.filter.bam"

    # Index
    samtools index "${BAM.getBaseName()}.filter.bam"
    mv "${BAM.getBaseName()}.filter.bam.bai" "${BAM.getBaseName()}.filter.bai"
    """
}
