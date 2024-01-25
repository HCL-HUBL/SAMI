process fastqc_trimmed {
	tag "${FASTQ.getSimpleName()}"

	cpus 1
	label 'monocore'
	label 'retriable'

	publishDir "${params.out}/QC/FastQC/trimmed", emit: "copy"

	when:
	!(FASTQ.name =~ /__input\.[0-9]+$/)

	input:
	path(FASTQ)

	output:
	path("${FASTQ.getSimpleName()}_fastqc.zip"), emit: QC_FASTQC_trimmed

	"""
	fastqc "$FASTQ" -o "."
	"""
}
