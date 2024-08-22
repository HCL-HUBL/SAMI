process rrna_interval {
	cpus 1
	time { 15.minute * task.attempt }
	memory { 500.MB * task.attempt }
	storeDir params.store

	input:
	path(GTF)
	path(chrNameLength)

	output:
	path("${GTF}.rRNA"), emit: rRNA

	"""
	# Header (consider unsorted to be safe)
	echo -e "@HD\tVN:1.0\tSO:unsorted" > "./${GTF}.rRNA"

	# Chromosomes (from STAR genome)
	sed -r 's/^(.+)\t(.+)\$/@SQ\tSN:\\1\tLN:\\2/' "${chrNameLength}" >> "./${GTF}.rRNA"

	# BED-like content
	grep -E 'transcript_(bio)?type "rRNA"' "$GTF" | awk -F "\t" '\$3 == "transcript" { id=gensub(/^.+transcript_id \"([^\"]+)\";.+\$/, "\\\\1", "g", \$9); print \$1"\t"\$4"\t"\$5"\t"\$7"\t"id }' >> "./${GTF}.rRNA"
	"""
}

process refflat {
	cpus 1
	time { 15.minute * task.attempt }
	memory { 4.GB * task.attempt }
	storeDir params.store

	input:
	path(GTF)

	output:
	path("${GTF}.refFlat"), emit: refFlat

	"""
	Rscript --vanilla "${projectDir}/scripts/gtfToRefFlat.R" "$GTF" "${GTF}.refFlat"
	"""
}

process rnaseqmetrics {
	tag "$sample"

	cpus 1
	time = { 1.hour * task.attempt }
	memory = { 4.GB * task.attempt }

	input:
	tuple val(sample), val(typeReads), path(BAM), path(BAI)
	val(typeGTF)
	path(refFlat)
	path(rRNA)
	val(strand)

	output:
	path("${sample}_${refFlat.name}_${typeGTF}.RNA_Metrics"), emit: RNA_Metrics

	"""
	# Run CollectRnaSeqMetrics
	java -Djava.io.tmpdir="${TMPDIR-/tmp/}" -Xmx4G -Duser.country=US -Duser.language=en -jar "\$picard" CollectRnaSeqMetrics \
		--INPUT $BAM \
		--OUTPUT "./${sample}_${refFlat.name}_${typeGTF}.RNA_Metrics" \
		--REF_FLAT "$refFlat" \
		--RIBOSOMAL_INTERVALS "$rRNA" \
		--STRAND_SPECIFICITY "${strand}" \
		--ASSUME_SORTED true
	"""
}
