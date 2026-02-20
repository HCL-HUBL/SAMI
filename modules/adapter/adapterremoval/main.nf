process adapterremoval {
	tag "$pair"

	cpus 2
	time { 10.minute * task.attempt }
	memory { 2.GB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(pair), val(type)

	output:
	path("${sample}_adapterremoval.log"), emit: log

	"""
	AdapterRemoval --identify-adapters --threads ${task.cpus} --file1 "$R1" --file2 "$R2" > "${sample}_adapterremoval.log"
	"""
}

process retrieveadapter {
	tag "$pair"

	cpus 1
	time { 5.minute * task.attempt }
	memory { 2.GB * task.attempt }
	publishDir "${params.output}/AdapterRemoval", mode: params.publish
	
	input:
	path('*')

	output:
	env(seqR1), emit: R1
	env(seqR2), emit: R2
	path("adapter.txt"), emit: adapter

	"""
	### Parse the log files from AdapterRemoval to get the sequences
	seqR1=$(awk '$0~/--adapter1:/ {print \$NF}' *_adapterremoval.log | uniq)
	seqR2=$(awk '$0~/--adapter2:/ {print \$NF}' *_adapterremoval.log | uniq)

	### Verify that only one adapter is present for each file
	if [ $(echo \$seqR1 | wc -w) -gt 1 ];  then echo "More than one adapter have been found for R1 files. Exit."; exit 1
	elif [ $(echo \$seqR2 | wc -w) -gt 1 ]; then echo "More than one adapter have been found for R2 files. Exit."; exit 1
	fi

	### Generate the file containing the adapter
	grep -h "\-\-adapter" *_adapterremoval.log | sort -u | sed -E 's/ +//' > adapter.txt
	"""
}
