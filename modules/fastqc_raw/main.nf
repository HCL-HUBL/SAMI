process fastqc_raw {

	cpus 1
	label 'monocore'
	label 'retriable'

	publishDir "${params.out}/QC/FastQC/raw", mode: "copy"

	when:
	!(FASTQ.name =~ /__input\.[0-9]+$/)

	input:
	path(FASTQ)

	output:
	path("${FASTQ.getSimpleName()}_fastqc.zip"), emit: QC_FASTQC_raw

	"""
	fastqc "$FASTQ" -o "."
	"""
}
