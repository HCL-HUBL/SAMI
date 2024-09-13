process splicing_depth {
	cpus { chunks }
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	path("BAM/*")
	path("BAM/*")
	path(BED)
	val(chunks)
	val(qMap)
	
	output:
	path("merged/depth.bed"), emit: BED

	"""
	# Split BED in X chunks
	mkdir in
	split -n l/${chunks} "${BED}" in/

	# Run samtools depth in parallel
	mkdir out
	declare -i i=0;
	for target in in/*
	do
	   samtools depth -a -H -Q ${qMap} -b "\$target" BAM/*.bam > out/\$(basename \$target) &
	   pids[\$i]=\$!
	   i=\$i+1
	done
	
	# Wait for processes to finish
	for pid in \${pids[*]}; do wait \$pid; done
	
	# Merge outputs
	mkdir merged
	head -1 -q out/* | uniq > "merged/depth.bed"
	grep -vh '#' out/* >> "merged/depth.bed"
	"""
}
