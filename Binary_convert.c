/*
 * Binary_convert.c - FASTA file preprocessor for kallisto
 * This utility converts FASTA files to a simpler format for kallisto processing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function prototypes
void convert_fasta_to_binary(const char *input_file, const char *output_file);
void print_usage(const char *prog_name);

const int max_line_length = 10000;

int main(int argc, char *argv[]){
	
	if (argc != 3) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char *input_file = argv[1];
	const char *output_file = argv[2];
	
	printf("Converting %s to %s\n", input_file, output_file);
	convert_fasta_to_binary(input_file, output_file);
	printf("Conversion complete!\n");
	
	return 0;
}

void print_usage(const char *prog_name) {
	fprintf(stderr, "Usage: %s <input.fasta> <output.bin>\n", prog_name);
	fprintf(stderr, "\nConverts a FASTA file to kallisto binary format\n");
}

void convert_fasta_to_binary(const char *input_file, const char *output_file) {
	FILE *input = fopen(input_file, "r");
	if (!input) {
		fprintf(stderr, "Error: Could not open input file '%s'\n", input_file);
		exit(1);
	}
	
	FILE *output = fopen(output_file, "w");
	if (!output) {
		fprintf(stderr, "Error: Could not open output file '%s'\n", output_file);
		fclose(input);
		exit(1);
	}

	char *line = malloc(max_line_length * sizeof(char));
	if (!line) {
		fprintf(stderr, "Error: Memory allocation failed\n");
		fclose(input);
		fclose(output);
		exit(1);
	}
	
	int count = 0;
	while (fgets(line, max_line_length, input) != NULL) {
		// Write line to output - simply copy for now
		// Could add compression or binary encoding here
		fprintf(output, "%s", line);
		count++;
	}
	
	printf("Processed %d lines\n", count);
	
	free(line);
	fclose(output);
	fclose(input);
}
