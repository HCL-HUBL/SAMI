process mutect2 {
    tag "$sample"

    cpus { params.CPU_mutect }
    label 'multicore'
    label 'nonRetriable'
    publishDir "${params.out}/variant", mode: params.publish
    scratch { params.scratch }

    input:
    tuple path(genomeFASTA), path(genomeFASTA_dict), path(genomeFASTA_fai)
    tuple path(gnomAD), path(gnomAD_index)
    tuple val(sample), val(type), path(BAM), path(BAI)

    output:
    tuple path("${sample}.filtered.vcf.gz"), path("${sample}.filtered.vcf.gz.tbi"), emit: filtered_VCF
    tuple path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), emit: unfiltered_VCF
    path("${sample}.unfiltered.vcf.gz.stats"), emit: stats

    """
    # Genomic subset
    if [ "${params.window}" = "" ]; then interval=""
    else                                 interval="--intervals ${params.window}"
    fi

    # Call variants
    gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" Mutect2 \
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
