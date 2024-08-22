process star_reindex {
	cpus 2
	time { 2.hour * task.attempt }
	memory { 40.GB + 10.GB * task.attempt }

	input:
	path(SJ)
	path(rawGenome)
	path(genomeGTF)
	path(dummy_R1)
	path(dummy_R2)
	val(genome)
	val(title)

	output:
	path("${genome}_${title}"), emit: genome

	"""
	mkdir -p "./reindex"
	STAR \
		--runThreadN 2 \
		--genomeDir "$rawGenome" \
		--readFilesIn "${dummy_R1}" "${dummy_R2}" \
		--sjdbFileChrStartEnd $SJ \
		--limitSjdbInsertNsj 5000000 \
		--sjdbInsertSave All \
		--sjdbGTFfile "$genomeGTF" \
		--outFileNamePrefix "./reindex/" \
		--outSAMtype None
	mv "./reindex/Log.out" "./reindex/_STARgenome/"
	mv "./reindex/_STARgenome/" "./${genome}_${title}/"
	"""
}
