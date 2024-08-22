// Collect and emit R1 + R2 pairs individually
def sample_sheet(sampleSheetPath) {
	FASTQ_list = []
	sampleSheet = file(sampleSheetPath)
	lines = sampleSheet.splitCsv(header: true)
	n = [:]
	for(line in lines) {
		// Pair ID
		if(n.containsKey(line["sample"])) {
			n[ line["sample"] ]++
		} else {
			n[ line["sample"] ] = 1
		}
		pair = line["sample"] + "_" + n[ line["sample"] ]
		
		// R1
		if(line["R1"] == "") {
			error "ERROR: A R1 file path is empty"
		} else {
			R1 = file(line["R1"])
		}
		
		// R2
		if(line["R2"] == "") {
			// Empty deterministic file (for -resume)
			R2_path = "${TMPDIR}/" + pair + "_empty-R2"
			R2 = file(R2_path)
			
			// Create file only if missing
			if(!R2.exists()) {
				R2.text = ''
			}
			
			type = "single"
		} else {
			R2 = file(line["R2"])
			type = "paired"
		}
		
		// Reformat as a list
		FASTQ_list << [ R1, R2, line["sample"], pair, type ]
	}
	
	// Reformat as a channel
	FASTQ = Channel.fromList(FASTQ_list)
	
	return FASTQ
}
