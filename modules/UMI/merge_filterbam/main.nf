process merge_filterbam {
    tag "$sample"

    cpus 2
    time { 1.hour * task.attempt }
	memory { 20.GB + 5.GB * task.attempt }

    input:
    tuple val(sample), val(type), path(BAM_mapped), path(BAM_unmapped), path(BAM_forUnmappedRead)
    tuple path(genomeFASTA), path(genomeFASTAdict), path(genomeFASTAindex)

    output:
    tuple val(sample), val(type), path("out/${sample}.DNA.bam"), emit: BAM

    """
    ### fgbio command
    fgBioExe="java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4g -jar \$fgbio"

    ### Sort the BAM by query name for gatk - VW need to be put before in STAR_pass2
    mkdir tmp
    java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" SortSam \
        -INPUT "${BAM_mapped}" \
        -OUTPUT mapped_sorted.bam \
        -SORT_ORDER queryname \
        --TMP_DIR "\$(pwd)/tmp"
    java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" SortSam \
        -INPUT "${BAM_unmapped}" \
        -OUTPUT unmapped_sorted.bam \
        -SORT_ORDER queryname \
        --TMP_DIR "\$(pwd)/tmp"

    ### Values from ROCHE pipeline
    gatk MergeBamAlignment \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_RETAIN RX \
        --ALIGNED_BAM mapped_sorted.bam \
        --UNMAPPED_BAM unmapped_sorted.bam \
        --OUTPUT "${sample}.temp.bam" \
        --REFERENCE_SEQUENCE "${genomeFASTA}" \
        --SORT_ORDER coordinate \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --CLIP_OVERLAPPING_READS false \
        --TMP_DIR "\$(pwd)/tmp"

    ### Values from Thomas (HCL_nextflow)
    samtools sort -n -u "${sample}.temp.bam" |
        \${fgBioExe} --async-io --compression 1 FilterConsensusReads \
        --ref=${genomeFASTA} \
        --input=/dev/stdin \
        --output="${sample}.filterConsensus.bam" \
        --min-reads=1 \
        --max-read-error-rate=0.05 \
        --max-base-error-rate=0.1 \
        --min-base-quality=1 \
        --max-no-call-fraction=0.2

    ### Get unmapped read from STAR_pass1 and put them after (from: https://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/)
    ### 0x4 4  UNMAP        0x108 264 MUNMAP,SECONDARY
    ### 0x8 8  MUNMAP       0x104 260 UNMAP,SECONDARY
    ### 0xc 12 UNMAP,MUNMAP 0x100 256 SECONDARY
    mkdir out
	samtools view -b -f4 -F264  "${BAM_forUnmappedRead}" > tmps1.bam
    samtools view -b -f8 -F260  "${BAM_forUnmappedRead}" > tmps2.bam
    samtools view -b -f12 -F256 "${BAM_forUnmappedRead}" > tmps3.bam
    samtools merge -o "out/${sample}.DNA.bam" "${sample}.filterConsensus.bam" "tmps"?".bam"
    """
}
