process splicing_annotation {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 4.GB * task.attempt }
	storeDir params.store

	input:
	path(genomeGTF)
	val(species)
	val(genome)
	val(chromosomes)
	
	output:
	path("${genomeGTF}.introns.rds"), emit: introns
	path("${genomeGTF}.exons.rdt"), emit: exons
	path("${genomeGTF}.genes.rdt"), emit: genes

	shell:
	template 'annotation.R'
}
