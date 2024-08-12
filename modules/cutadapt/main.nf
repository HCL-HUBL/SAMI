process cutadapt {
    tag "$sample"

	cpus { params.CPU_cutadapt }
	label 'monocore'
	label 'retriable'

	input:
	tuple path(R1), path(R2), val(sample), val(type)
	val(trimR1)
	val(trimR2)

	output:
	tuple path("${R1.getName().replaceFirst(/^(.+)([_\.]f(ast)?q(.gz)?)$/, '$1_cutadapt$2')}"), path("${R2.getName().replaceFirst(/^(.+)([_\.]f(ast)?q(.gz)?)$/, '$1_cutadapt$2')}"), val(sample), val(type), emit: FASTQ
	path("${sample}_cutadapt.log"), emit: log

	"""
	# Base cutadapt command
	command="cutadapt -j ${params.CPU_cutadapt} --minimum-length 20"
	
	# R1 adapter
	command="\${command} -a \"${trimR1}\""
	
	# R2 adapter
	if [[ ${type} == "paired" ]] && [[ ${trimR2} != "" ]]
	then
		command="\${command} -A \"${trimR2}\""
	fi
	
	# R1 output file
	R1_out="${R1.getName().replaceFirst(/^(.+)([_\.]f(ast)?q(.gz)?)$/, '$1_cutadapt$2')}"
	command="\${command} -o \"\$R1_out\""
	
	# R2 output file
	R2_out="${R2.getName().replaceFirst(/^(.+)([_\.]f(ast)?q(.gz)?)$/, '$1_cutadapt$2')}"
	if [[ ${type} == "paired" ]]
	then
		command="\${command} -p \"\$R2_out\""
	fi
	
	# R1 input file
	command="\${command} \"${R1}\""
	
	# R2 input file
	if [[ ${type} == "paired" ]]; then command="\${command} \"${R2}\""; fi
	
	# Execute
	\$command > "${sample}_cutadapt.log"
	"""
}
