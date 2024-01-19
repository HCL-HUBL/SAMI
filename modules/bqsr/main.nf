process BQSR {

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    publishDir "${params.out}/BQSR", mode: "copy"
    scratch { params.scratch }

    input:
    tuple(file(genomeFASTA), file(genomeFASTA_dict), file(genomeFASTA_fai))
    tuple(file(gnomAD), file(gnomAD_index))
    tuple(file(COSMIC), file(COSMIC_index))
    tuple(val(sample), val(type), file(BAM), file(BAI))

    output:
    tuple(val(sample), val(type), file("${BAM.getBaseName()}.BQSR.bam"), file("${BAM.getBaseName()}.BQSR.bai")), emit: BAM_BQSR
    file("${BAM.getBaseName()}.BQSR.clean"), emit: BQSR_clean

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
