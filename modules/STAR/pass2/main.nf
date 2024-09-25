process star_pass2 {
	tag "$sample"

	cpus 5
	time { 1.hour * task.attempt }
	memory { 30.GB + 5.GB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(type), val(RG)
	path(reindexedGenome)
	path(genomeGTF)
	val(protrude)
	val(multimap)
	
	output:
	tuple val(sample), val(type), path("${sample}.DNA.bam"), emit: BAM_DNA
	tuple val(sample), val(type), path("${sample}.RNA.bam"), emit: BAM_RNA
	tuple val(sample), val(type), path("${sample}.isize.txt"), emit: isize
	path("${sample}_Chimeric.out.junction"), emit: chimeric
	path("${sample}_Log.final.out"), emit: log

	"""
	# FASTQ files
	if [ "$type" = "paired" ];   then readFilesIn="\\"${R1.join(",")}\\" \\"${R2.join(",")}\\""
	elif [ "$type" = "single" ]; then readFilesIn="\\"${R1.join(",")}\\""
	else                         echo "Unknow type '$type'"; exit 1
	fi

	# Align
	mkdir -p "./$sample"
	STAR \
		--runThreadN ${task.cpus} \
		--twopassMode None \
		--genomeDir "$reindexedGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn \$readFilesIn \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "./${sample}/" \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--chimOutType Junctions WithinBAM \
		--quantMode TranscriptomeSAM \
		--outSAMattrRGline $RG \
		--sjdbGTFfile "$genomeGTF" \
		--alignEndsProtrude ${protrude} ConcordantPair \
		--alignInsertionFlush Right \
		--alignSJDBoverhangMin 4 \
		--alignSJstitchMismatchNmax 3 -1 3 3 \
		--alignSplicedMateMapLmin 16 \
		--alignSplicedMateMapLminOverLmate 0 \
		--chimJunctionOverhangMin 8 \
		--chimScoreJunctionNonGTAG -4 \
		--chimSegmentMin 10 \
		--chimMultimapNmax 1 \
		--chimNonchimScoreDropMin 10 \
		--outFilterMultimapNmax ${multimap} \
		--outFilterMismatchNmax 5 \
		--outSJfilterOverhangMin 16 8 8 8 \
		--outSJfilterDistToOtherSJmin 0 0 0 0 \
		--peOverlapNbasesMin 12 \
		--peOverlapMMp 0.1

	mv "./${sample}/Log.final.out" "./${sample}_Log.final.out"
	mv "./${sample}/SJ.out.tab" "./${sample}_SJ.out.tab"
	mv "./${sample}/Chimeric.out.junction" "./${sample}_Chimeric.out.junction"
	mv "./${sample}/Aligned.out.bam" "./${sample}.DNA.bam"
	mv "./${sample}/Aligned.toTranscriptome.out.bam" "./${sample}.RNA.bam"

	# Export ISIZE sample (empty in single-end)
	samtools view -f 0x2 -f 0x80 "./${sample}.RNA.bam" | cut -f9 | head -1000000 > "./${sample}.isize.txt"
	"""
}
