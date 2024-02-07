process splitn {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'nonRetriable'
    scratch { params.scratch }

    input:
    tuple path(genomeFASTA), path(genomeFASTA_dict), path(genomeFASTA_fai)
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    tuple val(sample), val(type), path("${BAM.getBaseName()}.splitN.bam"), path("${BAM.getBaseName()}.splitN.bai"), emit: BAM_splitN
    path("${BAM.getBaseName()}.splitN.clean"), emit: splitN_clean

    """
    gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" SplitNCigarReads \
        --input "$BAM" \
        --reference "$genomeFASTA" \
        --output "${BAM.getBaseName()}.splitN.bam" \
        --tmp-dir "."

    # Link input BAM for cleaning
    ln -s "${BAM.toRealPath()}" "${BAM.getBaseName()}.splitN.clean"
    """
}
