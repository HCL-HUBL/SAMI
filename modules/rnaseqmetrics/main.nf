process rnaseqmetrics {

    cpus 1
    label 'monocore'
    label 'retriable'
    publishDir "${params.out}/QC/rnaSeqMetrics", mode: "copy"

    input:
    tuple(val(sample), val(type), path(BAM), path(BAI), path(rRNA), path(refFlat))

    output:
    path("${sample}_${refFlat.name}_*.RNA_Metrics"), emit: QC_rnaSeqMetrics

    when:
    rRNA.name.replaceFirst(/\.rRNA$/, '') == refFlat.name.replaceFirst(/\.refFlat$/, '')

    """
    # GTF type
    if [[ "${refFlat.name.replaceFirst(/\.refFlat$/, '')}" == "${genomeGTF.getVal().name}" ]]
    then
    type="genome"
    elif [[ "${refFlat.name.replaceFirst(/\.refFlat$/, '')}" == "${targetGTF.getVal().name}" ]]
    then
    type="target"
    else
        echo "Unrecognized refFlat file"
    exit 1
    fi

    # Run CollectRnaSeqMetrics
    java -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CollectRnaSeqMetrics \
        INPUT=$BAM \
        OUTPUT="./${sample}_${refFlat.name}_${type}.RNA_Metrics" \
        REF_FLAT="$refFlat" \
        RIBOSOMAL_INTERVALS="$rRNA" \
        STRAND_SPECIFICITY="${params.stranded_Picard}" \
        ASSUME_SORTED=true
    """
}
