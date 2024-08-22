process indexfasta {
	cpus 1
	time { 5.minute * task.attempt }
	memory { 5.GB * task.attempt }
	storeDir params.store

	input:
	path(genomeFASTA)

	output:
	tuple path(genomeFASTA), path("${genomeFASTA.getBaseName()}.dict"), path("${genomeFASTA}.fai"), emit: indexedFASTA

	"""
	# Dictionnary
	java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CreateSequenceDictionary \
		--REFERENCE "$genomeFASTA" \
		--OUTPUT "${genomeFASTA.getBaseName()}.dict"
	# Index
	samtools faidx "$genomeFASTA"
	"""
}
