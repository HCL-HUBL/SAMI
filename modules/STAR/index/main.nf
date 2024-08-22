process star_index {
	cpus 10
	time { 2.hour * task.attempt }
	memory { 40.GB + 10.GB * task.attempt }
	storeDir params.store

	input:
	path(genomeFASTA)
	path(genomeGTF)

	output:
	path("${params.genome}_raw"), emit: genome
	path("${params.genome}_raw/chrNameLength.txt"), emit: chrom

	"""
	mkdir -p "./${params.genome}_raw"
	STAR \
		--runThreadN ${cpus} \
		--runMode genomeGenerate \
		--genomeDir "./${params.genome}_raw" \
		--genomeFastaFiles "$genomeFASTA" \
		--sjdbGTFfile "$genomeGTF"
	"""
}
