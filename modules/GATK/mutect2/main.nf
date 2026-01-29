process mutect2 {
	tag "$sample"

	cpus 4
	time { 1.hour * task.attempt }
	memory { 5.GB  * task.attempt }
	publishDir "${params.output}/variants", mode: params.publish

	input:
	tuple path(genomeFASTA), path(genomeFASTA_dict), path(genomeFASTA_fai)
	tuple path(germline), path(germline_index)
	tuple path(PoN), path(PoN_index)
	tuple val(sample), val(type), path(BAM), path(BAI)

	output:
	tuple val(sample), path("${sample}.filtered.vcf.gz"), path("${sample}.filtered.vcf.gz.tbi"), emit: filtered_VCF
	tuple val(sample), path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), emit: unfiltered_VCF
	path("${sample}.unfiltered.vcf.gz.stats"), emit: stats

	"""
	# Call variants
	gatk --java-options "-Xmx4G -Duser.country=US -Duser.language=en" Mutect2 \
		--input "$BAM" \
		--reference "$genomeFASTA" \
		--output "${sample}.unfiltered.vcf.gz" \
		--germline-resource "$germline" \
		--panel-of-normals "$PoN" \
		--native-pair-hmm-threads ${task.cpus}

	# Filter variants
	gatk FilterMutectCalls \
		--variant "${sample}.unfiltered.vcf.gz" \
		--reference "$genomeFASTA" \
		--output "${sample}.filtered.vcf.gz"
	"""
}
