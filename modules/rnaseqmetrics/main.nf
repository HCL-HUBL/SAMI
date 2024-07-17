process rnaseqmetrics {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'retriable'

    input:
    tuple val(sample), val(type), path(BAM), path(BAI), path(rRNA), path(refFlat)
    path(genomeGTF)
    path(targetGTF)

    output:
    path("${sample}_${refFlat.name}_*.RNA_Metrics"), emit: QC_rnaSeqMetrics

    when:
    rRNA.name.replaceFirst(/\.rRNA$/, '') == refFlat.name.replaceFirst(/\.refFlat$/, '')

    """
    # GTF type
    if [[ "${refFlat.name.replaceFirst(/\.refFlat$/, '')}" == "${genomeGTF.name}" ]]
    then
        typeGtf="genome"
    elif [[ "${refFlat.name.replaceFirst(/\.refFlat$/, '')}" == "${targetGTF.name}" ]]
    then
        typeGtf="target"
    else
        echo "Unrecognized refFlat file"
        exit 1
    fi

    # Run CollectRnaSeqMetrics
    java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CollectRnaSeqMetrics \
        --INPUT $BAM \
        --OUTPUT "./${sample}_${refFlat.name}_\${typeGtf}.RNA_Metrics" \
        --REF_FLAT "$refFlat" \
        --RIBOSOMAL_INTERVALS "$rRNA" \
        --STRAND_SPECIFICITY "${params.stranded_Picard}" \
        --ASSUME_SORTED true
    """
}
