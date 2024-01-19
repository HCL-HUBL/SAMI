process mutect2 {

    cpus { params.CPU_mutect }
    label 'multicore'
    label 'nonRetriable'
    publishDir "${params.out}/Mutect2", mode: "copy"
    scratch { params.scratch }

    input:
    tuple(file(genomeFASTA), file(genomeFASTA_dict), file(genomeFASTA_fai))
    tuple(file(gnomAD), file(gnomAD_index))
    tuple(val(sample), val(type), file(BAM), file(BAI))

    output:
    tuple(file("${sample}.filtered.vcf.gz"), file("${sample}.filtered.vcf.gz.tbi")), emit: filtered_VCF
    tuple(file("${sample}.unfiltered.vcf.gz"), file("${sample}.unfiltered.vcf.gz.tbi")), emit: unfiltered_VCF
    path("${sample}.unfiltered.vcf.gz.stats"), emit: Mutect2_stats

    """
    # Genomic subset
    if [ "${params.window}" = "" ]; then interval=""
    else                                 interval="--intervals ${params.window}"
    fi

    # Call variants
    gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" Mutect2 \$interval \
        --input "$BAM" \
        --reference "$genomeFASTA" \
        --output "${sample}.unfiltered.vcf.gz" \
        --germline-resource "$gnomAD" \
        --native-pair-hmm-threads ${params.CPU_mutect} \
        --independent-mates

    # Filter variants
    gatk FilterMutectCalls \
        --variant "${sample}.unfiltered.vcf.gz" \
        --reference "$genomeFASTA" \
        --output "${sample}.filtered.vcf.gz"
    """
}
