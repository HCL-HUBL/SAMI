#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
	
	int lineNumber = 1;    // Number of current STDIN line
	int lineIndex;         // Character position in current line
	int wordIndex;         // Character position in current word/cell
	int colIndex;          // Column number in current line
	char lineBuffer[32];   // To store the whole current line
	char wordBuffer[32];   // To store the current word/cell
	char chrom[16];        // Chromosome of current line
	char prevChrom[16];    // Chromosome of previous line
	int pos;               // Genomic position of current line
	int prevPos;           // Genomic position of previous line
	int rleStart;          // Genomic start of current RLE
	int depth;             // Sequencing depth of current line
	int prevDepth;         // Sequencing depth of previous line
	
	while(fgets(lineBuffer, sizeof(lineBuffer), stdin) != NULL) {
		lineIndex = 0;
		wordIndex = 0;
		colIndex = 0;

		// Parse current line
		while(lineBuffer[lineIndex] != '\0') {
//			printf("line %i, character %i = [%c]\n", lineNumber, lineIndex, lineBuffer[lineIndex]);
			if(lineBuffer[lineIndex] == '\t' || lineBuffer[lineIndex] == '\n') {
				// End of word
				wordBuffer[wordIndex] = '\0';
				
//				printf("line %i, column %i, word=[%s]\n", lineNumber, colIndex, wordBuffer);
				
				// Interpret the word
				if(colIndex == 0)        {
					if(wordIndex > sizeof(wordBuffer)) { fprintf(stderr, "Chromosome column too large line %i\n", lineNumber);
					} else                             { strncpy(chrom, wordBuffer, sizeof(chrom));
					}
				} else if(colIndex == 1) { pos = atoi(wordBuffer);
				} else if(colIndex == 2) { depth = atoi(wordBuffer);
				} else                   { fprintf(stderr, "Too many columns line %i\n", lineNumber);
				}
				
				// Reset word index
				wordIndex = 0;
					
				// Next column
				colIndex++;
			} else {
				// Bufferize current character
				wordBuffer[wordIndex] = lineBuffer[lineIndex];
				
				// Next character in word
				wordIndex++;
			}
			
			// Next character in line
			lineIndex++;
		}

//		printf("line %i, chrom=[%s], pos=[%i], depth=[%i]\n", lineNumber, chrom, pos, depth);
		
		if(lineNumber == 1) {
			// First line starts a new RLE
			rleStart = pos;
		} else {
			// Process current line
			if(strcmp(chrom, prevChrom) == 0 && depth == prevDepth && pos == prevPos + 1) {
				// Same RLE as previous line
			} else {
				// Output previous RLE (chrom / start / end / depth)
				printf("%s\t%i\t%i\t%i\n", prevChrom, rleStart, prevPos, prevDepth);
				
				// Start a new RLE
				rleStart = pos;
			}
		}
	
		// Remember current line as previous line
		strncpy(prevChrom, chrom, sizeof(prevChrom));
		prevPos = pos;
		prevDepth = depth;
		
		lineNumber++;
	}
	
	// Output last RLE
	printf("%s\t%i\t%i\t%i\n", prevChrom, rleStart, prevPos, prevDepth);
	
	return EXIT_SUCCESS;
}

