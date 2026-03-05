process cutadapt {
	tag "$pair"

	cpus 2
	time { 20.minute * task.attempt }
	memory { 2.GB * task.attempt }

	input:
	tuple path(R1), path(R2), val(sample), val(pair), val(type)
	val(trimR1)
	val(trimR2)

	output:
	tuple path("out/${R1.getName()}"), path("out/${R2.getName()}"), val(sample), val(pair), val(type), emit: FASTQ
	path("${pair}_cutadapt.log"), emit: log

	"""
	# Base cutadapt command
	command="cutadapt -j ${task.cpus} --minimum-length 20"

	# R1 adapter
	command="\${command} -a \"${trimR1}\""

	# R2 adapter
	if [[ ${type} == "paired" ]] && [[ ${trimR2} != "" ]]
	then
		command="\${command} -A \"${trimR2}\""
	fi

	# R1 output file
	mkdir -p "out"
	R1_out="out/${R1.getName()}"
	command="\${command} -o \"\$R1_out\""

	# R2 output file
	mkdir -p "out"
	R2_out="out/${R2.getName()}"
	if [[ ${type} == "paired" ]]
	then
		command="\${command} -p \"\$R2_out\""
	else
		touch "\$R2_out"
	fi

	# R1 input file
	command="\${command} \"${R1}\""

	# R2 input file
	if [[ ${type} == "paired" ]]; then command="\${command} \"${R2}\""; fi

	# Execute
	\$command > "${pair}_cutadapt.log"
	"""
}
