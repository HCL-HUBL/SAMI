// Collect and emit R1 + R2 pairs individually
def sample_sheet(sampleSheetPath) {
	FASTQ_list = []
	sampleSheet = file(sampleSheetPath)
	lines = sampleSheet.splitCsv(header: true)
	n = [:]
	for(line in lines) {
		// R1
		if(line["R1"] == "") {
			error "ERROR: A R1 file path is empty"
		} else {
			R1 = file(line["R1"])
		}
		
		// R2
		if(line["R2"] == "") {
			// Empty deterministic file (for -resume)
			if(n.containsKey(line["sample"])) {
				n[ line["sample"] ]++
			} else {
				n[ line["sample"] ] = 1
			}
			R2_path = "${TMPDIR}/" + line["sample"] + "_empty-R2_" + n[ line["sample"] ]
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
		FASTQ_list << [ R1, R2, line["sample"], type ]
	}
	
	// Reformat as a channel
	FASTQ = Channel.fromList(FASTQ_list)
	
	return FASTQ
}

/*
def sample_sheet(sampleSheetPath) {
    // Collect R1 and R2 per sample
	FASTQ_map = [:];
	sampleSheet = file(sampleSheetPath)
	lines = sampleSheet.splitCsv(header: true)
	for(line in lines) {
		sampleName = line["sample"]
		if(FASTQ_map.containsKey(sampleName)) {
			FASTQ_map[sampleName]["R1"] << line["R1"]
			FASTQ_map[sampleName]["R2"] << line["R2"]
		} else {
			FASTQ_map[sampleName] = [ "R1": [ line["R1"] ], "R2": [ line["R2"] ] ]
		}
	}
	
	// Check content
	FASTQ_list = [];
	FASTQ_map.each { sampleName, sample ->
		// Check lists
		if(sample["R1"].size() != sample["R2"].size()) {
			error "ERROR: R1 and R2 file counts differ for $sampleName"
		}
		if(sample["R1"].size() == 0) {
			error "ERROR: No R1 file for $sampleName"
		}
		
		// Check elements
		anyPE = false
		for(int i = 0; i < sample["R1"].size(); i++) {
			if(sample["R1"][i] == "") {
				error "ERROR: Empty R1 file path for $sampleName"
			} else {
				sample["R1"][i] = file(sample["R1"][i])
			}
			if(sample["R2"][i] == "") {
				sample["R2"][i] = file("mktemp".execute().text.replaceAll("\\s",""))
				if(anyPE) {
					error "ERROR: Mixed single and paired end files for $sampleName"
				}
			} else {
				anyPE = true
				sample["R2"][i] = file(sample["R2"][i])
			}
		}
		
		// Paired-end or single-end
		if(anyPE) { type = "paired"
		} else    { type = "single"
		}
		
		// Reformat as a list
		FASTQ_list << [ "R1": sample["R1"], "R2": sample["R2"], "sample": sampleName, "type": type ]
	}
	
	// Reformat as a channel
	FASTQ = Channel.fromList(FASTQ_list)
	
	return FASTQ
}
*/