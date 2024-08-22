process umi_consensus{
	tag "$sample"

	cpus { params.CPU_umi }
	time { 1.hour * task.attempt }
	memory { 5.GB * task.attempt }

	input:
	tuple val(sample), path(BAM), val(type), val(RG)

	output:
	tuple val(sample), path("${sample}_family_size_histogram.txt"), emit: histogram
	tuple path("${sample}.consensus_R1.fastq.gz"), path("${sample}.consensus_R2.fastq.gz"), val(sample), val(type), env(newRG), emit: FASTQ
	tuple val(sample), path("${sample}.consensus.bam"), emit: BAM_unmapped

	"""
	set -eo pipefail

	### Create a new RG at sample level
	newRG=\$(echo "$RG" | sed -E 's/ *, *.+\$//')
	newRG_ID="consensus"
	newRG=\$(echo "\$newRG" | sed -E "s/ID:([^\t]+)/ID:\$newRG_ID/")

	### fgbio command
	fgBioExe="java -Djava.io.tmpdir="\${TMPDIR-/tmp/}" -Xmx4g -XX:-UsePerfData -jar \$fgbio"

	\${fgBioExe} --async-io --compression 0 CopyUmiFromReadName \
		--input="${BAM}" \
		--output=/dev/stdout \
		| \${fgBioExe} --async-io --compression 0 SortBam \
		--input=/dev/stdin \
		--output=/dev/stdout \
		--sort-order=Queryname \
		| \${fgBioExe} --async-io --compression 0 SetMateInformation \
		--input=/dev/stdin \
		--output=/dev/stdout \
		| \${fgBioExe} --async-io --compression 0 GroupReadsByUmi \
		--input=/dev/stdin \
		--output=/dev/stdout \
		--family-size-histogram="${sample}_family_size_histogram.txt" \
		--raw-tag=RX \
		--assign-tag=MI \
		--strategy=Adjacency \
		--edits=1 \
		| \${fgBioExe} --async-io CallMolecularConsensusReads \
		--input=/dev/stdin \
		--output="${sample}.consensus.bam" \
		--error-rate-pre-umi 30 \
		--error-rate-post-umi 30 \
		--output-per-base-tags true \
		--min-reads 1 \
		--max-reads 50 \
		--min-input-base-quality 10 \
		--read-name-prefix="csr" \
		--threads "${params.CPU_umi}" \
		--read-group-id="\${newRG_ID}"

	### Convert into FASTQ
	if [ "$type" = "paired" ]
	then
		# Paired-end
		samtools collate -u -O "${sample}.consensus.bam" | \
			samtools fastq \
			-1 "${sample}.consensus_R1.fastq.gz" \
			-2 "${sample}.consensus_R2.fastq.gz" \
			-n
	else
		# Single-end
		samtools collate -u -O "${sample}.consensus.bam" | \
			samtools fastq \
			-0 "${sample}.consensus_R1.fastq.gz" \
			-n
		touch "${sample}.consensus_R2.fastq.gz"
	fi
	"""
}
