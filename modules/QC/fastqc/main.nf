process fastqc {
	tag "${FASTQ.getSimpleName()}"

	cpus 1
	label 'monocore'
	label 'retriable'

	input:
	path(FASTQ)

	output:
	path("*_fastqc.zip"), emit: zip

	"""
	fastqc "$FASTQ" -o "."
	"""
}
