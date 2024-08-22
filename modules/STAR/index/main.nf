process star_index {
	cpus 10
	time { 2.hour * task.attempt }
	memory { 40.GB + 10.GB * task.attempt }
	storeDir params.store

	input:
	path(genomeFASTA)
	path(genomeGTF)
	val(genome)

	output:
	path("${genome}_raw"), emit: genome
	path("${genome}_raw/chrNameLength.txt"), emit: chrom

	"""
	mkdir -p "./${genome}_raw"
	STAR \
		--runThreadN ${task.cpus} \
		--runMode genomeGenerate \
		--genomeDir "./${genome}_raw" \
		--genomeFastaFiles "$genomeFASTA" \
		--sjdbGTFfile "$genomeGTF"
	"""
}
