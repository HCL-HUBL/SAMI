process bcftools {
	tag "$sample"

	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/variants", mode: params.publish

	input:
	tuple val(sample), path(VCF), path(TBI)
	val format
	val include
	val exclude

	output:
	tuple val(sample), path("${VCF}.tsv"), emit: TSV

	"""
	if [ "$format" == "full" ]
	then
		# Extract field names from VCF header
		FORMAT="\$(zcat "$VCF" | grep -E '^##FORMAT=<ID=([^,]+),.+\$' | sed -E 's/^##FORMAT=<ID=([^,]+),.+\$/[%\\1]/' | tr '\n' '\t')"
		INFO="\$(zcat "$VCF" | grep -E '^##INFO=<ID=([^,]+),.+\$' | sed -E 's/^##INFO=<ID=([^,]+),.+\$/%\\1/' | tr '\n' '\t')"
		
		# Add fixed fields
		final_format="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t\${FORMAT}\${INFO::-1}\n"
	else
		final_format="$format"
	fi
	
	# Apply filtering
	if [ -z "$include" ]
	then
		if [ -z "$exclude" ]
		then
			bcftools query -HH --format "\$final_format" "$VCF" > "${VCF}.tsv"
		else
			bcftools query -HH --format "\$final_format" --exclude "$exclude" "$VCF" > "${VCF}.tsv"
		fi
	else
		if [ -z "$exclude" ]
		then
			bcftools query -HH --format "\$final_format" --include "$include" "$VCF" > "${VCF}.tsv"
		else
			bcftools query -HH --format "\$final_format" --include "$include" --exclude "$exclude" "$VCF" > "${VCF}.tsv"
		fi
	fi
	"""
}
