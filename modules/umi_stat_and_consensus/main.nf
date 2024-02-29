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
    ### Create a temporary directory for fgbio execution
    tmpdir="${sample}_\$(date +%H%M%S)"
    mkdir "\${tmpdir}"

    ### fgbio command
    fgBioExe="java -Xmx4g -Djava.io.tmpdir=\${tmpdir} -jar \$fgbio"

    ### Function to run at the end to clean the temporary files
    function cleanup()
    {
        rm -rf "${sample}.changeName.bam" "${sample}.copy.bam" "${sample}.sort.bam" "${sample}.mate.bam" "${sample}.grpUmi.bam" "\${tmpdir}"
    }

    ### Clean the temporary file when the program exit
    trap cleanup EXIT

    ### Get only the ID of the RG
    newRG=\$(echo "${RG}" | cut -f1 | sed 's/ID://')

    ### Change the "_" into a ":" before the UMI in read name
    samtools view -h "${BAM}" | sed -r 's/(^[^\t]*:[0-9]*)_([ATCGN]*)\t/\\1:\\2\t/' | samtools view -b > "${sample}.changeName.bam"

    ### Put UMI as a tag in the bam file
    \${fgBioExe} --async-io --compression 1 CopyUmiFromReadName \
        --input="${sample}.changeName.bam" \
        --output="${sample}.copy.bam"

    ### Put mate info after sorting
    \${fgBioExe} --async-io --compression 1 SortBam \
        --input="${sample}.copy.bam" \
        --output="${sample}.sort.bam" \
        --sort-order=Queryname
    \${fgBioExe} --async-io --compression 1 SetMateInformation \
        --input="${sample}.sort.bam" \
        --output="${sample}.mate.bam"

    ### Group reads per UMI
    \${fgBioExe} --async-io --compression 1 GroupReadsByUmi \
        --input="${sample}.mate.bam" \
        --output="${sample}.grpUmi.bam" \
        --family-size-histogram="${sample}_family_size_histogram.txt" \
        --raw-tag=RX \
        --assign-tag=MI \
        --strategy=Adjacency \
        --edits=1

    ### Get consensus
    \${fgBioExe} --async-io CallMolecularConsensusReads \
        --input="${sample}.grpUmi.bam" \
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
