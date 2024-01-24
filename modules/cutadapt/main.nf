process cutadapt {

	cpus { params.CPU_cutadapt }
	label 'monocore'
	label 'retriable'

	publishDir "${params.out}/cutadapt", mode: "copy"

	input:
	tuple(path(R1), path(R2), val(sample), val(type), val(RG))

	output:
	tuple(path(env(outR1)), path(env(outR2)), val(sample), val(type), val(RG)), emit: FASTQ_STAR1
	path(env(outR1)), emit: R1_trimmed
	path(env(outR2)), emit: R2_trimmed
	path(env(outQC)), emit: QC_cutadapt

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

	outR1="${R1.getSimpleName()}_cutadapt.fastq.gz"
	outR2="${R2.getSimpleName()}_cutadapt.fastq.gz"
	outQC="${sample}_cutadapt.log"
	"""
}
