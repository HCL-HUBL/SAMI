#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_BUFFER 1024

// Process the minima of bin sizes from the single argument
int process_args(
		char* raw_sizes,   // Minima of all bin sizes separated by commas, eg: '0,30,60,90,120,150'
		int** min_size,    // Where to store the parsed values
		int* n_sizes       // Where to store the amount of parsed values
	) {
	
	// At most half of string size elements
	int maxElements = strlen(raw_sizes) / 2;
	int* int_sizes = (int*) calloc(maxElements, sizeof(int));
	
	// Split and cast provided sizes
	int n = 0;
	char* token = strtok(raw_sizes, ",");
	while(token != NULL) {
		// Process one value
		int_sizes[n] = atoi(token);
		token = strtok(NULL, ",");
		
		// Check values
		if(int_sizes[n] < 0) {
			fprintf(stderr, "Expected positive sizes (%i)\n", int_sizes[n]);
			return EXIT_FAILURE;
		}
		if(n > 0 && int_sizes[n] <= int_sizes[n-1]) {
			fprintf(stderr, "Expected strictly increasing sizes (%i,%i)\n", int_sizes[n], int_sizes[n-1]);
			return EXIT_FAILURE;
		}
		
		n++;
	}
	
	// Export
	*n_sizes = n;
	*min_size = int_sizes;
	
	return EXIT_SUCCESS;
}

// Handles the header lines of both FASTQ files
int process_header(
		char* R1_line,     // Line buffer for R1
		char* R2_line,     // Line buffer for R2
		char* R1_header,   // Where to copy R1_line
		char* R2_header,   // Where to copy R2_line
		int line           // Number of currently processed FASTQ line
	) {
	
	if(R1_line[0] != '@') {
		fprintf(stderr, "Expected a '@' character in R1 file line %i\n", line);
		return EXIT_FAILURE;
	}
	if(R2_line[0] != '@') {
		fprintf(stderr, "Expected a '@' character in R2 file line %i\n", line);
		return EXIT_FAILURE;
	}
	strncpy(R1_header, R1_line, LINE_BUFFER);
	strncpy(R2_header, R2_line, LINE_BUFFER);
	
	return EXIT_SUCCESS;
}

// Handles the sequence lines of both FASTQ files
int process_sequence(
		char* R1_line,    // Line buffer for R1
		char* R2_line,    // Line buffer for R2
		char* R1_seq,     // Where to copy R1_line
		char* R2_seq,     // Where to copy R2_line
		int* R1_length,   // Where to store the length of R1_line (with new line characters)
		int* R2_length    // Where to store the length of R2_line (with new line characters)
	) {
	
	*R1_length = strlen(R1_line);
	*R2_length = strlen(R2_line);
	strncpy(R1_seq, R1_line, LINE_BUFFER);
	strncpy(R2_seq, R2_line, LINE_BUFFER);
	
	return EXIT_SUCCESS;
}

// Handles the "+" lines of both FASTQ files
int process_separator(
		char* R1_line,   // Line buffer for R1
		char* R2_line,   // Line buffer for R1
		int line         // Number of currently processed FASTQ line
	) {
	
	if(R1_line[0] != '+' || strlen(R1_line) != 2) {
		fprintf(stderr, "Expected a single (%i) '+' character in R1 file line %i\n", strlen(R1_line), line);
		return EXIT_FAILURE;
	}
	if(R2_line[0] != '+' || strlen(R1_line) != 2) {
		fprintf(stderr, "Expected a single '+' character in R2 file line %i\n", line);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}

// Handles the quality lines of both FASTQ files
int process_quality(
		char* R1_line,               // Line buffer for R1
		char* R2_line,               // Line buffer for R2
		int R1_length,               // Where to store the length of R1_line (with new line characters)
		int R2_length,               // Where to store the length of R2_line (with new line characters)
		char* R1_header,             // Header of current R1 previously encountered
		char* R2_header,             // Header of current R2 previously encountered
		char* R1_seq,                // Sequence of current R1 previously encountered
		char* R2_seq,                // Sequence of current R2 previously encountered
		FILE** R1_output_handles,    // File handles of output R1 FASTQ files
		FILE** R2_output_handles,    // File handles of output R2 FASTQ files
		FILE* R1_discarded_handle,   // File handle of R1 discarded reads (R1 and R2 lengths differ too much)
		FILE* R2_discarded_handle,   // File handle of R2 discarded reads (R1 and R2 lengths differ too much)
		int n_sizes,                 // Amount of bins
		int* min_size,               // Minimal value of each bin
		int* exported,               // Amount of exported read pairs for each bin
		int* discarded_diff,         // Where to increment the amount of discarded read pairs (R1 and R2 lengths differ too much)
		int* reads,                  // Where to increment the amount of processed reads so far
		int line                     // Number of currently processed FASTQ line
	) {
	
	int diff_length;
	int min_length;
	
	// Check consistency
	if(R1_length != strlen(R1_line)) {
		fprintf(stderr, "Size inconsistency between sequence and quality in R1 file line %i\n", line);
		return EXIT_FAILURE;
	}
	if(R2_length != strlen(R2_line)) {
		fprintf(stderr, "Size inconsistency between sequence and quality in R2 file line %i\n", line);
		return EXIT_FAILURE;
	}
	
	// True length (without line breaks)
	R1_length = strcspn(R1_line, "\r\n");
	R2_length = strcspn(R2_line, "\r\n");
	
	// Read length
	if(R1_length > R2_length) {
		diff_length = R1_length - R2_length;
		min_length = R2_length;
	} else {
		diff_length = R2_length - R1_length;
		min_length = R1_length;
	}
	
	FILE* R1_handle = NULL;
	FILE* R2_handle = NULL;
	if(diff_length <= 5) {
		for(int i = n_sizes-1; i >= 0; i--) {
			if(min_length >= min_size[i]) {
				// Handle to export to
				R1_handle = R1_output_handles[i];
				R2_handle = R2_output_handles[i];
				
				// Successfully exported
				exported[i]++;
				break;
			}
		}
	} else {
		// Handle to export to
		R1_handle = R1_discarded_handle;
		R2_handle = R2_discarded_handle;
		
		(*discarded_diff)++;
	}
	
	if(R1_handle == NULL || R2_handle == NULL) {
		fprintf(stderr, "Read size (%i) out of range line %i\n", min_length, line);
		return EXIT_FAILURE;
	}
	
	// Write R1
	fputs(R1_header, R1_handle);
	fputs(R1_seq, R1_handle);
	fputs("+\n", R1_handle);
	fputs(R1_line, R1_handle);
	
	// Write R2
	fputs(R2_header, R2_handle);
	fputs(R2_seq, R2_handle);
	fputs("+\n", R2_handle);
	fputs(R2_line, R2_handle);
	
	(*reads)++;
	
	return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
	
	if(argc != 4) {
		fprintf(stderr, "USAGE : BinReads R1.fastq R2.fastq MIN_SIZE_1[,MIN_SIZE_2[...]]\n", argc);
		return EXIT_FAILURE;
	}
	
	char R1_buffer[LINE_BUFFER];   // To store one line from R1 FASTQ
	char R2_buffer[LINE_BUFFER];   // To store one line from R2 FASTQ
	char* R1_line;                 // R1_buffer or NULL
	char* R2_line;                 // R2_buffer or NULL
	char R1_header[LINE_BUFFER];   // To store the header line of R1
	char R2_header[LINE_BUFFER];   // To store the header line of R2
	char R1_seq[LINE_BUFFER];      // To store the sequence line of R1
	char R2_seq[LINE_BUFFER];      // To store the sequence line of R2
	int line = 0;                  // Line number inside the R1 and R2 files
	int reads = 0;                 // Reads fully processed so far
	int readLine = 0;              // Line number inside a read (1 = header, 2 = sequence, 3 = '+', 4 = quality)
	int R1_length;                 // Size of R1 (including new line characters)
	int R2_length;                 // Size of R2 (including new line characters)
	int status;                    // Exit codes returned by functions
	int discarded_diff = 0;        // Amount of discarded read pairs (R1 and R2 lengths differ too much)
	
	// Process arguments
	int* min_size;     // Minimal value of each bin
	int n_sizes = 0;   // Amount of bins
	status = process_args(
		argv[3],
		&min_size,
		&n_sizes
	);
	if(status == EXIT_FAILURE) return EXIT_FAILURE;
	
	// Amount of exported read pairs for each bin
	int* exported = (int*) calloc(n_sizes, sizeof(int));
	
	// Open input R1 file handler
	FILE* R1_input_handle = fopen(argv[1], "r");
	if(R1_input_handle == NULL) {
		fprintf(stderr, "Error opening file \"%s\"\n", argv[1]);
		return EXIT_FAILURE;
	}
	
	// Open input R2 file handler
	FILE* R2_input_handle = fopen(argv[2], "r");
	if(R2_input_handle == NULL) {
		fprintf(stderr, "Error opening file \"%s\"\n", argv[2]);
		return EXIT_FAILURE;
	}
	
	// Open output file handlers
	char file_name[1024];   // To store dynamic file names to open
	FILE** R1_output_handles = calloc(n_sizes, sizeof(FILE*));
	FILE** R2_output_handles = calloc(n_sizes, sizeof(FILE*));
	for(int i=0; i<n_sizes; i++) {
		// R1 output
		sprintf(file_name, "R1_gte%i.fastq", min_size[i]);
		R1_output_handles[i] = fopen(file_name, "w");
		if(R1_output_handles[i] == NULL) {
			fprintf(stderr, "Error creating file \"%s\"\n", file_name);
			return EXIT_FAILURE;
		}
		
		// R2 output
		sprintf(file_name, "R2_gte%i.fastq", min_size[i]);
		R2_output_handles[i] = fopen(file_name, "w");
		if(R2_output_handles[i] == NULL) {
			fprintf(stderr, "Error creating file \"%s\"\n", file_name);
			return EXIT_FAILURE;
		}
	}
	
	// R1 discarded output
	sprintf(file_name, "R1_discarded.fastq");
	FILE* R1_discarded_handle = fopen(file_name, "w");
	if(R1_discarded_handle == NULL) {
		fprintf(stderr, "Error creating file \"%s\"\n", file_name);
		return EXIT_FAILURE;
	}
	
	// R2 discarded output
	sprintf(file_name, "R2_discarded.fastq");
	FILE* R2_discarded_handle = fopen(file_name, "w");
	if(R2_discarded_handle == NULL) {
		fprintf(stderr, "Error creating file \"%s\"\n", file_name);
		return EXIT_FAILURE;
	}
	
	// Loop over FASTQ file lines
	while(1) {
		
		line++;
		
		// Get one line from each FASTQ
		R1_line = fgets(R1_buffer, sizeof(R1_buffer), R1_input_handle);
		R2_line = fgets(R2_buffer, sizeof(R1_buffer), R2_input_handle);
		
		// Break loop
		if(R1_line) {
			if(R2_line) {
				// Everything is fine, continue
			} else {
				fprintf(stderr, "R1 file seems longer than R2 file\n");
				return EXIT_FAILURE;
			}
		} else {
			if(R2_line) {
				fprintf(stderr, "R2 file seems longer than R1 file\n");
				return EXIT_FAILURE;
			} else {
				break;
			}
		}
		
		// Move inside the read paragraph
		readLine++;
		if(readLine == 5) readLine = 1;
		
		// Check FASTQ format
		switch(readLine) {
			
			case 1:
				// Header line
				status = process_header(
					R1_line,
					R2_line,
					R1_header,
					R2_header,
					line
				);
				if(status == EXIT_FAILURE) return EXIT_FAILURE;
				break;
			
			case 2:
				// Sequence line
				status = process_sequence(
					R1_line,
					R2_line,
					R1_seq,
					R2_seq,
					&R1_length,
					&R2_length
				);
				if(status == EXIT_FAILURE) return EXIT_FAILURE;
				break;
			
			case 3:
				// Separator line
				status = process_separator(
					R1_line,
					R2_line,
					line
				);
				if(status == EXIT_FAILURE) return EXIT_FAILURE;
				break;
			
			case 4:
				// Quality line
				status = process_quality(
					R1_line,
					R2_line,
					R1_length,
					R2_length,
					R1_header,
					R2_header,
					R1_seq,
					R2_seq,
					R1_output_handles,
					R2_output_handles,
					R1_discarded_handle,
					R2_discarded_handle,
					n_sizes,
					min_size,
					exported,
					&discarded_diff,
					&reads,
					line
				);
				if(status == EXIT_FAILURE) return EXIT_FAILURE;
				
				break;
		}
	}
	
	// Close files
	fclose(R1_input_handle);
	fclose(R2_input_handle);
	for(int i = 0; i < n_sizes; i++) {
		fclose(R1_output_handles[i]);
		fclose(R2_output_handles[i]);
	}
	fclose(R1_discarded_handle);
	fclose(R2_discarded_handle);
	
	// Summary
	fprintf(stdout, "Pairs processed = %i\n", reads);
	fprintf(stdout, "Pairs exported :\n");
	for(int i = 0; i < n_sizes; i++) fprintf(stdout, "- >= %i bp = %i\n", min_size[i], exported[i]);
	fprintf(stdout, "Pairs discarded (diff > 5) = %i\n", discarded_diff);
	
	return EXIT_SUCCESS;
}

