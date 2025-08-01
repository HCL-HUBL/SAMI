process uncompress {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 2.GB * task.attempt }

	input:
	path(file)
	
	output:
	path("out/*")

	"""
	shopt -s nocasematch
	
	mkdir out
	filename="$file"
	if [[ "\$filename" =~ \\.gz\$ ]]
	then
		# Uncompress
		zcat "\$filename" > "out/\${filename%.gz}"
	else
		# Leave unchanged
		ln -s "\$filename" out/
	fi
	"""
}
