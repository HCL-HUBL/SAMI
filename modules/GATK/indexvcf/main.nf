process indexvcf {
	cpus 1
	time { 1.hour * task.attempt }
	memory { 2.GB  * task.attempt }
	storeDir params.store

	input:
	path(gnomAD)
	path(PoN)

	output:
	tuple path("$gnomAD"), path("${gnomAD}.tbi"), emit: germline
	tuple path("$PoN"), path("${PoN}.tbi"), emit: PoN

	"""
	tabix "$gnomAD"
	tabix "$PoN"
	"""
}
