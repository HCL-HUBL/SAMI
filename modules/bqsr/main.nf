process bqsr {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/BQSR", mode: "copy"
    scratch { params.scratch }

    input:
    tuple path(genomeFASTA), path(genomeFASTA_dict), path(genomeFASTA_fai)
    tuple path(gnomAD), path(gnomAD_index)
    tuple path(COSMIC), path(COSMIC_index)
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    tuple val(sample), val(type), path("${BAM.getBaseName()}.BQSR.bam"), path("${BAM.getBaseName()}.BQSR.bai"), emit: BAM_BQSR
    path("${BAM.getBaseName()}.BQSR.clean"), emit: BQSR_clean

    """
    # Compute model
    gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" BaseRecalibrator \
        --input "$BAM" \
        --reference "$genomeFASTA" \
        --known-sites "$gnomAD" \
        --known-sites "$COSMIC" \
        --output "${sample}.BQSR" \
        --tmp-dir "."

    # Apply model
    gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" ApplyBQSR \
        --input "$BAM" \
        --reference "$genomeFASTA" \
        --bqsr-recal-file "${sample}.BQSR" \
        --output "${BAM.getBaseName()}.BQSR.bam" \
        --tmp-dir "."

    # Link input BAM for cleaning
    ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.BQSR.clean"
    """
}
