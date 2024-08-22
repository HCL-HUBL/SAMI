process annotation {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 2.GB * task.attempt }
	storeDir params.store

	input:
	path(genomeGTF)

	output:
	path("${genomeGTF}.introns.rds"), emit: introns
	path("${genomeGTF}.exons.rdt"), emit: exons
	path("${genomeGTF}.genes.rdt"), emit: genes

	"""
	Rscript --vanilla "${projectDir}/scripts/annotation.R" "$genomeGTF" "$params.species" "$params.genome" "$params.chromosomes"
	"""
}
