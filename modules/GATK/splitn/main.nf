process splitn {
    tag "$sample"

    cpus 1
	time { 1.hour * task.attempt }
	memory { 5.GB * task.attempt }

    input:
    tuple path(genomeFASTA), path(genomeFASTA_dict), path(genomeFASTA_fai)
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    tuple val(sample), val(type), path("${BAM.getBaseName()}.splitN.bam"), path("${BAM.getBaseName()}.splitN.bai"), emit: BAM

    """
    gatk --java-options "-Djava.io.tmpdir=\"\${TMPDIR-/tmp/}\" -Xmx4G -Duser.country=US -Duser.language=en" SplitNCigarReads \
        --input "$BAM" \
        --reference "$genomeFASTA" \
        --output "${BAM.getBaseName()}.splitN.bam" \
        --tmp-dir "."
    """
}
