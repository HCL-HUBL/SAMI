process bam_sort {
	tag "$sample"

	cpus 4
	time { 1.hour * task.attempt }
	memory { 4.GB * task.attempt }
	publishDir "${params.output}/BAM", mode: params.publish

	input:
	tuple val(sample), val(type), path(BAM)

	output:
	tuple val(sample), val(type), path("${BAM.getBaseName()}.sort.bam"), path("${BAM.getBaseName()}.sort.bai"), emit: BAM

	"""
	# Sort
	samtools sort -o "${BAM.getBaseName()}.sort.bam" -T ./${sample} -@ 3 "$BAM"

	# Index
	samtools index "${BAM.getBaseName()}.sort.bam"
	mv "${BAM.getBaseName()}.sort.bam.bai" "${BAM.getBaseName()}.sort.bai"
	"""
}
