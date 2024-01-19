process clean_dna_bam {

    cpus 1
    label 'nonRetriable'

    // Never scratch
    scratch false
    stageInMode 'symlink'
    executor 'local'

    when:
    params.clean_BAM

    input:
    path(clean)

    """
    # BAM file
    BAM="\$(readlink "$clean")"
    echo -n '' > \$BAM

    # BAI file
    BAI="\${BAM%.bam}.bai"
    if [ -f "\$BAI" ]; then echo -n '' > \$BAI; fi
    """
}
