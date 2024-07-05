process umi_stat_and_consensus{
    tag "$sample"

    cpus { params.CPU_umi }
    label 'multicore'
    label 'retriable'

    input:
    tuple val(sample), path(BAM), val(type), val(RG)
    tuple path(R1), path(R2), val(sample), val(type), val(RG)

    output:
    tuple val(sample), path("${sample}_family_size_histogram.txt"), emit: UMI_stat
    path("${sample}_family_size_histogram.txt"), emit: UMI_histo
    tuple path("${sample}.consensus_R1.fastq.gz"), path("${sample}.consensus_R2.fastq.gz"), val(sample), val(type), val(RG), emit: FASTQ_STAR2
    tuple val(sample), path("${sample}.consensus.bam"), emit: BAM_unmapped

    """
    set -eo pipefail

    ### fgbio command
    fgBioExe="java -Djava.io.tmpdir="\${TMPDIR-/tmp/}" -Xmx4g -jar \$fgbio"

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
        --read-group-id="\${newRG}" \
        --threads "${params.CPU_umi}"

    ### Convert into FASTQ
    samtools collate -u -O "${sample}.consensus.bam" | \
        samtools fastq \
        -1 "${sample}.consensus_R1.fastq.gz" \
        -2 "${sample}.consensus_R2.fastq.gz" \
        -n
    """
}
