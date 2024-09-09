#include <stdio.h>
#include "htslib/htslib/sam.h"

/* Chained list of gaps (on a specific chromosome)
   - hts_pos_t start             Start position of the gap (smaller coordinate)
   - hts_pos_t end               End position of the gap (bigger coordinate)
   - unsigned int reads          Count of reads which showed this gap
   - struct gap_list *previous   Pointer to next node in (sorted by coordinate) list [TODO: use it !]
   - struct gap_list *next       Pointer to previous node in (sorted by coordinate) list
*/
struct gap_list {
	hts_pos_t start;
	hts_pos_t end;
	unsigned int reads;
	struct gap_list *previous;
	struct gap_list *next;
};

/* Add an element to the chained list of gaps
   - struct gap_list *start_gap   Pointer to the first element of the list
   - struct gap_list *next_gap    Pointer to the list element suppoed to be after the new one
   - char beyond                  0 means 'next_gap' is after the new element, 1 means 'next_gap' is already the last element but insertion must be done after
   - hts_pos_t start              Start genomic position of the new element
   - hts_pos_t end                End genomic position of the new element
*/
struct gap_list * insert_gap(struct gap_list *start_gap, struct gap_list *next_gap, char beyond, hts_pos_t start, hts_pos_t end) {
	// Create new node
	struct gap_list *new = malloc(sizeof(struct gap_list));
	new -> start = start;
	new -> end = end;
	new -> reads = 1;
	
	if(start_gap == NULL) {
		// Initialize list
		new -> previous = NULL;
		new -> next = NULL;
		return new;
	} else if(beyond) {
		// New last element
		new -> previous = next_gap;
		new -> next = NULL;
		next_gap -> next = new;
		return start_gap;
	} else if(start_gap == next_gap) {
		// New first element
		new -> previous = NULL;
		new -> next = next_gap;
		next_gap -> previous = new;
		return new;
	} else {
		// New intersticial element
		new -> previous = next_gap -> previous;
		new -> next = next_gap;
		next_gap -> previous -> next = new;
		next_gap -> previous = new;
		return start_gap;
	}
}

int main(int argc, char *argv[])
{
	// Arguments
	if(argc != 2) {
		fprintf(stderr, "Expecting 1 argument\n");
		return 1;
	}
	
	// Minimum mapping quality for a read to be considered
	uint8_t min_qmap = 20;
	
	// Open BAM file
	bam1_t *b = bam_init1();
	samFile *in = sam_open(argv[1], "r");
	if (in == NULL) return -1;
	
	// Get SAM header
	sam_hdr_t *header = sam_hdr_read(in);
	if (header == NULL) return -1;
	
	// Array of gap lists (one per chromosome)
	struct gap_list **gaps = calloc(header -> n_targets, sizeof(struct gap_list *));
	
	// Loop over all reads
	int i = 0;
	while (sam_read1(in, header, b) >= 0) {
		// Progression
		i++;
		if(i % 100000 == 0) fprintf(stderr, "%d reads processed\n", i);
		
		// Current read "core" structure
		bam1_core_t *core = &b -> core;
		
		// Mapping quality filter
		if(core -> qual < min_qmap) {
			continue;
		}
		
		// Filter out secondary alignments
		if((core -> flag & BAM_FSECONDARY) != 0) {
			continue;
		}
		
		// Alignment start position (0-based to 1-based)
		hts_pos_t posBefore = core -> pos + 1;
		hts_pos_t posAfter;
		
		// Loop over CIGAR operations
		unsigned int *cigar = bam_get_cigar(b);
		for (int k = 0; k < core -> n_cigar; ++k) {
			// CIGAR operation code (char)
			char opchr = bam_cigar_opchr(cigar[k]);
			
			// CIGAR operation length (in bp)
			int oplen = bam_cigar_oplen(cigar[k]);
			
			// Does the CIGAR operation consume the reference ?
			int op = bam_cigar_op(cigar[k]);
			int optype = bam_cigar_type(op);
			// int qcons = optype & 1;
			int rcons = (optype & 2) / 2;
			
			// Genomic position at the end of the operation
			if(rcons) { posAfter = posBefore + oplen;
			} else    { posAfter = posBefore;
			}
			
			// Alignment gap
			if(opchr == 'N') {
//				printf("%c %d : %s:%ld-%ld\n", opchr, oplen, chr, posBefore, posAfter);
				
				struct gap_list *current_gap = gaps[ core -> tid ];
				
				if (current_gap == NULL) {
					// New node is the new list start
					gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], NULL, 0, posBefore, posAfter);
				} else {
					while(1) {
						if(current_gap -> start > posBefore) {
							// Insert new node before current node
							gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 0, posBefore, posAfter);
							break;
						} else if(current_gap -> start == posBefore) {
							if(current_gap -> end > posAfter) {
								// Insert new node before current node
								gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 0, posBefore, posAfter);
								break;
							} else if(current_gap -> end == posAfter) {
								// Exact match, increment
								current_gap -> reads ++;
								break;
							}
						}
						
						if (current_gap -> next == NULL) {
							// Insert new node at the end of the list, after current node
							gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 1, posBefore, posAfter);
							break;
						}
						
						// Move on
						current_gap = current_gap -> next;
					}
				}
			}
			
			// Move to prepare next operation
			posBefore = posAfter;
		}
	}
	
	// Print gap list as BED (end - 1 as in STAR's *.SJ.out.tab)
	for (int k = 0; k < header -> n_targets; ++k) {
		char* chrom = header -> target_name[k];
		struct gap_list *current_gap = gaps[k];
		while(current_gap != NULL) {
			printf("%s\t%ld\t%ld\t%d\n", chrom, current_gap -> start, current_gap -> end - 1, current_gap -> reads);
			current_gap = current_gap -> next;
		}
	}
	
	// Cleanup
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
}
