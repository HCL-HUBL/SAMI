process star_pass1 {
	tag "$sample"

	cpus 5
	time { 1.hour * task.attempt }
	memory { 30.GB + 5.GB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(type), val(RG)
	path(rawGenome)
	path(genomeGTF)
	val(protrude)

	output:
	path("${sample}.SJ.out.tab"), emit: junctions
	tuple val(sample), path("${sample}.pass1.bam"), val(type), val(RG), emit: BAM_DNA
	path("${sample}_Log.final.out"), emit: log

	"""
	mkdir -p "./$sample"

	# FASTQ files
	if [ "$type" = "paired" ];   then readFilesIn="\\"${R1.join(",")}\\" \\"${R2.join(",")}\\""
	elif [ "$type" = "single" ]; then readFilesIn="\\"${R1.join(",")}\\""
	else                         echo "Unknow type '$type'"; exit 1
	fi

	STAR \
		--runThreadN ${task.cpus} \
		--twopassMode None \
		--genomeDir "$rawGenome" \
		--genomeLoad NoSharedMemory \
		--readFilesIn \$readFilesIn \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "./" \
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
		--outFilterMultimapNmax 3 \
		--outFilterMismatchNmax 5 \
		--outSJfilterOverhangMin 16 8 8 8 \
		--outSJfilterDistToOtherSJmin 0 0 0 0 \
		--peOverlapNbasesMin 12 \
		--peOverlapMMp 0.1

	mv ./SJ.out.tab ./${sample}.SJ.out.tab
	mv "./Aligned.out.bam" "./${sample}.pass1.bam"
	mv "./Log.final.out" "./${sample}_Log.final.out"
	"""
}
