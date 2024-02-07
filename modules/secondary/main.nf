process secondary {
    tag "$sample"

    cpus 1
    label 'monocore'
    label 'retriable'

    input:
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    path("${sample}_mqc.yaml"), emit: QC_secondary

    """
    # Total read count (fast, from index)
    total=\$(samtools idxstats $BAM | awk 'BEGIN{x=0} {x=x+\$3+\$4} END{print x}')

    # Unaligned
    una=\$(samtools view -c -f 0x4 $BAM)

    # Secondary alignments
    sec=\$(samtools view -c -f 0x100 $BAM)

    # Primary alignments (deduced)
    pri=\$(bc <<< "scale=2; \$total-\$una-\$sec")

    # MultiQC regular file header
    cat <<-EOF > "./${sample}_mqc.yaml"
    id: 'Secondary_section'
    section_name: 'Secondary alignments'
    description: 'as a proportion of all alignments (reads) returned by STAR. A read aligning in multiple locations is duplicated in the BAM file (one entry for each alignment), the proportion of secondary alignments shows to which extent the read count was artificially increased by this phenomenon.'
    plot_type: 'bargraph'
    pconfig:
    id: 'Secondary_bargraph'
    title: 'Secondary alignments'
    data:
    ${sample}: {Primary: \${pri}, Unaligned: \${una}, Secondary: \${sec}}
    EOF
    """
}
