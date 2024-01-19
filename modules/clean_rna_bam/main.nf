process clean_rna_bam {

    cpus 1
    label 'nonRetriable'

    // Never scratch
    scratch false
    stageInMode 'symlink'
    executor 'local'

    when:
    params.clean_BAM && !params.RNA_BAM

    input:
    tuple(val(sample), val(type), path(BAM))

    """
    echo -n '' > "\$(readlink "$BAM")"
    """
}
