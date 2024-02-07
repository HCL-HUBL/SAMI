process cutadapt {
    tag "$sample"

	cpus { params.CPU_cutadapt }
	label 'monocore'
	label 'retriable'

	input:
	tuple path(R1), path(R2), val(sample), val(type), val(RG)

	output:
	tuple path("${R1.getSimpleName()}_cutadapt.fastq.gz"), path("${R2.getSimpleName()}_cutadapt.fastq.gz"), val(sample), val(type), val(RG), emit: FASTQ_STAR1
    path("${R1.getSimpleName()}_cutadapt.fastq.gz"), emit: R1_trimmed
	path("${R2.getSimpleName()}_cutadapt.fastq.gz"), emit: R2_trimmed
	path("${sample}_cutadapt.log"), emit: QC_cutadapt

	"""
	if [[ ${params.trimR1} != "" ]] && [[ ${params.trimR2} != "" ]]
	then
	cutadapt -j ${params.CPU_cutadapt} \
		-a "${params.trimR1}" \
		-A "${params.trimR2}" \
		--minimum-length 20 \
		-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
		-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
		"${R1}" "${R2}" > "${sample}_cutadapt.log"
	elif [[ ${params.trimR1} != "" ]]
	then
	cutadapt -j ${params.CPU_cutadapt} \
		-a "${params.trimR1}" \
		--minimum-length 20 \
		-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
		-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
		"${R1}" "${R2}" > "${sample}_cutadapt.log"
	else
		cutadapt -j ${params.CPU_cutadapt} \
		-A "${params.trimR2}" \
		--minimum-length 20 \
		-o "${R1.getSimpleName()}_cutadapt.fastq.gz" \
		-p "${R2.getSimpleName()}_cutadapt.fastq.gz" \
		"${R1}" "${R2}" > "${sample}_cutadapt.log"
	fi
	"""
}
