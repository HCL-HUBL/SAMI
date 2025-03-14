process filterduplicates {
	tag "$sample"

	cpus 1
	time { 1.hour * task.attempt }
	memory { 5.GB * task.attempt }

	input:
	tuple val(sample), val(type), path(BAM), path(BAI)

	output:
	tuple val(sample), val(type), path("${BAM.getBaseName()}.filter.bam"), path("${BAM.getBaseName()}.filter.bai"), emit: BAM

	"""
	# Filter
	samtools view -b -F 0x400 "$BAM" -o "${BAM.getBaseName()}.filter.bam"

	# Index
	samtools index "${BAM.getBaseName()}.filter.bam"
	mv "${BAM.getBaseName()}.filter.bam.bai" "${BAM.getBaseName()}.filter.bai"
	"""
}
