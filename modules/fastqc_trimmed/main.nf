process fastqc_trimmed {

	cpus 1
	label 'monocore'
	label 'retriable'

	publishDir "${params.out}/QC/FastQC/trimmed", emit: "copy"

	when:
	!(FASTQ.name =~ /__input\.[0-9]+$/)

	input:
	path(FASTQ)

	output:
	path(env(outQC)), emit QC_FASTQC_trimmed

	if(params.trimR1 != '' || params.trimR2 != '') {
		"""
		fastqc "$FASTQ" -o "."
		outQC="${FASTQ.getSimpleName()}_fastqc.zip"
		"""
	} else {
		// Bypass FastQC_trimmed
		"""
		outQC="${baseDir}/in/dummy.tsv"
		"""
	}
}
