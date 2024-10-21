process fastq_skip {
	tag "$sample"

	cpus 1
	time { 5.minute * task.attempt }
	memory { 500.MB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(pair), val(type)
	val(CN)
	val(PL)
	val(PM)

	output:
	tuple path(R1), path(R2), val(sample), val("${type.first()}"), stdout, emit: FASTQ

	shell:
	template 'fastq_skip.R'
}
