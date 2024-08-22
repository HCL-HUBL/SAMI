process fastqc {
	tag "${FASTQ.getName()}"

	cpus 1
	time { 1.hour * task.attempt }
	memory { 4.GB * task.attempt }

	input:
	path(FASTQ)

	output:
	path("*_fastqc.zip"), emit: zip

	"""
	fastqc "$FASTQ" -o "."
	"""
}
