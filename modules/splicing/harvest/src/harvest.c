#include <stdio.h>
#include "htslib/htslib/sam.h"

/* Pair orientations (0 for single-end)
   - F = Forward alignment
   - R = Reverse alignment
   - 1 = R1 (either on the left, genomically speaking, or the right)
   - 2 = R2 (either on the left, genomically speaking, or the right)
   - e.g. F2R1 means R2 is on the left aligning forward and R1 on the right aligning reverse (plausible orientation)
*/
#define F1R2 1
#define F2R1 2
#define F1F2 3
#define F2F1 4
#define R1R2 5
#define R2R1 6
#define R1F2 7
#define R2F1 8

/* Chained list of gaps (on a specific chromosome)
   - hts_pos_t start             Start position of the gap (smaller coordinate)
   - hts_pos_t end               End position of the gap (bigger coordinate)
   - hts_pos_t end               End position of the gap (bigger coordinate)
   - unsignedchar orientation    Pair orientation (using integer codes defined above)
   - struct gap_list *previous   Pointer to next node in (sorted by coordinate) list [TODO: use it !]
   - struct gap_list *next       Pointer to previous node in (sorted by coordinate) list
*/
struct gap_list {
	hts_pos_t start;
	hts_pos_t end;
	unsigned char orientation;
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
   - unsignedchar orientation     Pair orientation (using integer codes defined above)
*/
struct gap_list * insert_gap(struct gap_list *start_gap, struct gap_list *next_gap, char beyond, hts_pos_t start, hts_pos_t end, unsigned char orientation) {
	// Create new node
	struct gap_list *new = malloc(sizeof(struct gap_list));
	new -> start = start;
	new -> end = end;
	new -> orientation = orientation;
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

/* Determine pair orientation from a single read
   - bam1_core_t *read   Pointer to the read core structure
*/
unsigned char get_orientation(bam1_core_t *read) {
	// Is read R1 or R2 ?
	char isReadR1;
	if((read -> flag & BAM_FREAD1) != 0)        { isReadR1 = 1;
	} else if((read -> flag & BAM_FREAD2) != 0) { isReadR1 = 0;
	} else                                      { return 0;
	}
	
	// Is read first genomically ?
	if(read -> tid != read -> mtid) return 0;
	char isReadFirst = (read -> pos < read -> mpos);	
	
	// Read and mate strands
	char isReadReverse = ((read -> flag & BAM_FREVERSE) != 0);
	char isMateReverse = ((read -> flag & BAM_FMREVERSE) != 0);
	
	unsigned char output;
	if(isReadR1) {
		if(isReadFirst) {
			if(isReadReverse) {
				if(isMateReverse) { output = R1R2;
				} else            { output = R1F2;
				}
			} else {
				if(isMateReverse) { output = F1R2;
				} else            { output = F1F2;
				}
			}
		} else {
			if(isReadReverse) {
				if(isMateReverse) { output = R2R1;
				} else            { output = F2R1;
				}
			} else {
				if(isMateReverse) { output = R2F1;
				} else            { output = F2F1;
				}
			}
		}
	} else {
		if(isReadFirst) {
			if(isReadReverse) {
				if(isMateReverse) { output = R2R1;
				} else            { output = R2F1;
				}
			} else {
				if(isMateReverse) { output = F2R1;
				} else            { output = F2F1;
				}
			}
		} else {
			if(isReadReverse) {
				if(isMateReverse) { output = R1R2;
				} else            { output = F1R2;
				}
			} else {
				if(isMateReverse) { output = R1F2;
				} else            { output = F1F2;
				}
			}
		}
	}
	
	return output;
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
		
		// Pair orientation
		unsigned char orientation = get_orientation(core);
		
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
					gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], NULL, 0, posBefore, posAfter, orientation);
				} else {
					while(1) {
						if(current_gap -> start > posBefore) {
							// Insert new node before current node
							gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 0, posBefore, posAfter, orientation);
							break;
						} else if(current_gap -> start == posBefore) {
							if(current_gap -> end > posAfter) {
								// Insert new node before current node
								gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 0, posBefore, posAfter, orientation);
								break;
							} else if(current_gap -> end == posAfter) {
								if(current_gap -> orientation > orientation) {
									// Insert new node before current node
									gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 0, posBefore, posAfter, orientation);
									break;
								} else if(current_gap -> orientation == orientation) {
									// Exact match, increment
									current_gap -> reads ++;
									break;
								}
							}
						}
						
						if (current_gap -> next == NULL) {
							// Insert new node at the end of the list, after current node
							gaps[ core -> tid ] = insert_gap(gaps[ core -> tid ], current_gap, 1, posBefore, posAfter, orientation);
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
	
	// Orientation dictionnary
	char* dic[9];
	dic[0] = "unknown";
	dic[F1R2] = "F1R2";
	dic[F2R1] = "F2R1";
	dic[F1F2] = "F1F2";
	dic[F2F1] = "F2F1";
	dic[R1R2] = "R1R2";
	dic[R2R1] = "R2R1";
	dic[R1F2] = "R1F2";
	dic[R2F1] = "R2F1";
	
	// Print gap list as BED (end - 1 as in STAR's *.SJ.out.tab)
	for (int k = 0; k < header -> n_targets; ++k) {
		char* chrom = header -> target_name[k];
		struct gap_list *current_gap = gaps[k];
		while(current_gap != NULL) {
			printf(
				"%s\t%ld\t%ld\t%s\t%d\n",
				chrom,
				current_gap -> start,
				current_gap -> end - 1,
				dic[ current_gap -> orientation ],
				current_gap -> reads
			);
			current_gap = current_gap -> next;
		}
	}
	
	// Cleanup
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
}
