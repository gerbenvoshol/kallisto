/*
 * Kallisto C Implementation
 * A reimplementation of kallisto RNA-seq quantification tool in C
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <getopt.h>
#include <zlib.h>
#include <time.h>
#include "kseq.h"

// Declare kt_for() function from kthread.c
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

// Initialize kseq for FILE* reading
KSEQ_INIT(gzFile, gzread)

#define VERSION "0.3.0"

typedef struct{
	char * kmer_seq; // the kmer sequence
	int * t_labels; // list of labels of transcripts that contain the k-mer
	int from_reads;
} kmer_labeled;

typedef struct{
	int eq_class_id; // concatenation of all of the transcript names in the equivalency class
	int count; // the number of fragments(kmers) with the specified equivalency class
	int * eq_class_labels; // list of labels of transcripts that contain the k-mer
} eq_c_count;

// Function prototypes
uint32_t kmer_hash(char * kmer_seq);
uint32_t eq_hash(char * eq_class);
void store_name_idx(char * name, int t_name_idx);
char** get_kmers(char* seq, int k, int *size);
void create_htable(char * tscripts_file);
void store_read_counts(char * reads_file, int max_reads);
void store_paired_read_counts(char * reads_file1, char * reads_file2, int max_reads);
int * get_equiv_class(char* kmer);
int * compute_eq_class_intersection(int **eq_classes, int num_classes);
void increment_eq_class_count(int *eq_class);
void store_eqiv_classes();
uint32_t murmurhash (const char *key, uint32_t len, uint32_t seed);
char *clean(char *str);
void print_eqc_arr();
void print_t_lengths();
void free_hashTable();
void compare(double *probs);
void print_usage(const char *prog_name);
void print_version();
void write_abundances(const char *output_file, double *probs);
void write_abundances_with_bootstrap(const char *output_file, double *probs, double **bootstrap_results, int n_bootstrap);

//EM related
void EM(int *lengths, double eps);
void run_EM_iteration(int *lengths, double eps, double *probs_out);
void bootstrap_EM(int *lengths, double eps);
void update_parameters(double *alpha, double *probs, int *lengths, double **transcirpt_prob);
void update_alphas(double *alpha, double **transcirpt_prob);
void update_prob(double *alpha, double *probs, int *lengths);
void update_trans_probs(double *probs, double **transcirpt_prob);
double **alloc_matrix(int n, int r);
double likelihood(double *alpha, int *lengths);

// Data structure for parallel computation
typedef struct {
	double *alpha;
	double **transcript_prob;
	double *probs;
	int *lengths;
} em_data_t;

// Global variables
kmer_labeled **hashTable;
eq_c_count *eqc_arr;
int tr_name_chars = 100;  // Increased for longer names
int max_hasht_rows = 0;  // Will be set dynamically based on k-mer count
int max_per_row = 10;
int k = 31;  // Default k-mer size (kallisto default)
int num_t = 125;
int num_eq_classes;
int max_line_length = 10000;  // Increased for longer sequences
int num_eqc_found;
int bootstrap_samples = 0;  // Number of bootstrap samples (0 = no bootstrap)
int use_full_debruijn = 0;  // Use full de Bruijn graph (all k-mers) for pseudoalignment (0 = first k-mer only)

char ** t_names;
int * t_lengths;
char * g_output_file = NULL;  // Global output file path
int n_threads = 1;  // Default number of threads

void print_usage(const char *prog_name) {
	fprintf(stderr, "kallisto C implementation v%s\n\n", VERSION);
	fprintf(stderr, "Usage: %s [options] -i <index> -o <output>\n\n", prog_name);
	fprintf(stderr, "Required arguments:\n");
	fprintf(stderr, "  -i, --index <file>        Transcriptome index/reference file (FASTA/FASTQ format)\n");
	fprintf(stderr, "  -o, --output <file>       Output file for abundance estimates\n\n");
	fprintf(stderr, "Single-end mode:\n");
	fprintf(stderr, "  -r, --reads <file>        Input reads file (FASTA/FASTQ format, gzipped supported)\n\n");
	fprintf(stderr, "Paired-end mode:\n");
	fprintf(stderr, "  -1, --reads1 <file>       First reads file (FASTA/FASTQ format, gzipped supported)\n");
	fprintf(stderr, "  -2, --reads2 <file>       Second reads file (FASTA/FASTQ format, gzipped supported)\n\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "  -k, --kmer-size <int>     K-mer size (default: 31)\n");
	fprintf(stderr, "  -e, --epsilon <float>     EM convergence threshold (default: 0.01)\n");
	fprintf(stderr, "  -m, --max-reads <int>     Maximum number of reads to process (default: all)\n");
	fprintf(stderr, "  -t, --threads <int>       Number of threads (default: 1)\n");
	fprintf(stderr, "  -b, --bootstrap <int>     Number of bootstrap samples (default: 0)\n");
	fprintf(stderr, "  -d, --full-debruijn       Use full de Bruijn graph (all k-mers) for pseudoalignment\n");
	fprintf(stderr, "  -h, --help                Display this help message\n");
	fprintf(stderr, "  -v, --version             Display version information\n\n");
}

void print_version() {
	printf("kallisto C implementation version %s\n", VERSION);
}

int main(int argc, char *argv[])
{
	// Seed random number generator for bootstrap
	srand(time(NULL));
	
	char *index_file = NULL;
	char *reads_file = NULL;
	char *reads_file1 = NULL;
	char *reads_file2 = NULL;
	char *output_file = NULL;
	int max_reads = INT_MAX;
	double epsilon = 0.01;
	int paired_end = 0;
	
	static struct option long_options[] = {
		{"index",      required_argument, 0, 'i'},
		{"reads",      required_argument, 0, 'r'},
		{"reads1",     required_argument, 0, '1'},
		{"reads2",     required_argument, 0, '2'},
		{"output",     required_argument, 0, 'o'},
		{"kmer-size",  required_argument, 0, 'k'},
		{"epsilon",    required_argument, 0, 'e'},
		{"max-reads",  required_argument, 0, 'm'},
		{"threads",    required_argument, 0, 't'},
		{"bootstrap",  required_argument, 0, 'b'},
		{"full-debruijn", no_argument,    0, 'd'},
		{"help",       no_argument,       0, 'h'},
		{"version",    no_argument,       0, 'v'},
		{0, 0, 0, 0}
	};
	
	int opt;
	int option_index = 0;
	
	while ((opt = getopt_long(argc, argv, "i:r:1:2:o:k:e:m:t:b:dhv", long_options, &option_index)) != -1) {
		switch (opt) {
			case 'i':
				index_file = optarg;
				break;
			case 'r':
				reads_file = optarg;
				break;
			case '1':
				reads_file1 = optarg;
				paired_end = 1;
				break;
			case '2':
				reads_file2 = optarg;
				paired_end = 1;
				break;
			case 'o':
				output_file = optarg;
				break;
			case 'k':
				k = atoi(optarg);
				if (k <= 0 || k > 100) {
					fprintf(stderr, "Error: k-mer size must be between 1 and 100\n");
					return 1;
				}
				break;
			case 'e':
				epsilon = atof(optarg);
				if (epsilon <= 0) {
					fprintf(stderr, "Error: epsilon must be positive\n");
					return 1;
				}
				break;
			case 'm':
				max_reads = atoi(optarg);
				if (max_reads <= 0) {
					fprintf(stderr, "Error: max-reads must be positive\n");
					return 1;
				}
				break;
			case 't':
				n_threads = atoi(optarg);
				if (n_threads <= 0) {
					fprintf(stderr, "Error: threads must be positive\n");
					return 1;
				}
				break;
			case 'b':
				bootstrap_samples = atoi(optarg);
				if (bootstrap_samples < 0) {
					fprintf(stderr, "Error: bootstrap samples must be non-negative\n");
					return 1;
				}
				break;
			case 'd':
				use_full_debruijn = 1;
				break;
			case 'h':
				print_usage(argv[0]);
				return 0;
			case 'v':
				print_version();
				return 0;
			default:
				print_usage(argv[0]);
				return 1;
		}
	}
	
	// Check required arguments
	if (!index_file || !output_file) {
		fprintf(stderr, "Error: Missing required arguments\n\n");
		print_usage(argv[0]);
		return 1;
	}
	
	// Check read file arguments
	if (paired_end && (!reads_file1 || !reads_file2)) {
		fprintf(stderr, "Error: Both -1 and -2 must be specified for paired-end mode\n\n");
		print_usage(argv[0]);
		return 1;
	}
	
	if (!paired_end && !reads_file) {
		fprintf(stderr, "Error: Either -r or both -1 and -2 must be specified\n\n");
		print_usage(argv[0]);
		return 1;
	}
	
	if (paired_end && reads_file) {
		fprintf(stderr, "Error: Cannot use both -r and -1/-2 simultaneously\n\n");
		print_usage(argv[0]);
		return 1;
	}
	
	// Set global output file
	g_output_file = output_file;
	
	printf("Kallisto C implementation v%s\n", VERSION);
	printf("Index file:  %s\n", index_file);
	if (paired_end) {
		printf("Reads file 1: %s\n", reads_file1);
		printf("Reads file 2: %s\n", reads_file2);
		printf("Mode: Paired-end\n");
	} else {
		printf("Reads file:  %s\n", reads_file);
		printf("Mode: Single-end\n");
	}
	printf("Output file: %s\n", output_file);
	printf("K-mer size:  %d\n", k);
	printf("Epsilon:     %f\n", epsilon);
	printf("Threads:     %d\n", n_threads);
	if (bootstrap_samples > 0) {
		printf("Bootstrap:   %d samples\n", bootstrap_samples);
	}
	printf("Pseudoalignment: %s\n", use_full_debruijn ? "Full de Bruijn graph (all k-mers)" : "First k-mer only");
	printf("\n");

	create_htable(index_file);

	if (paired_end) {
		store_paired_read_counts(reads_file1, reads_file2, max_reads);
	} else {
		store_read_counts(reads_file, max_reads);
	}
	store_eqiv_classes();

	// Freeing the hash table because is not needed for the EM
	free_hashTable();

	printf("Equivalence classes: \n");
	print_eqc_arr();

	printf("\nRunning EM algorithm...\n");
	EM(t_lengths, epsilon);

	return 0;
}


// returns the index in the hash table that the k-mer sequence and associated transctipt labels are stored at
// use hash function from https://github.com/jwerle/murmurhash.c
// kallisto uses the same murmurhash3 function
uint32_t kmer_hash(char * kmer_seq)
{
	uint32_t seed = 5;
	uint32_t h_val =  murmurhash(kmer_seq, (uint32_t) strlen(kmer_seq), seed) % max_hasht_rows;

	return h_val;
}

// This method returns a list of the possible kmers of a transcripts
// *seq is the transcript, k is the length of the kmer and size if a
// variable used to know the amout of kmers in the transcript
char** get_kmers(char* seq, int k, int *size)
{
	int seq_len = strlen(seq);
	if (seq_len <= k) {
		*size = 0;
		return NULL;
	}
	int num_kmers = seq_len - k + 1;  // Fixed: should be +1, not just -k
	char ** kmers = malloc(num_kmers * sizeof(char*));
	if (!kmers) {
		fprintf(stderr, "Failed to allocate memory for kmers\n");
		*size = 0;
		return NULL;
	}

	for (int i = 0; i < num_kmers; i++){
		kmers[i] = malloc((k + 1) * sizeof(char));  // +1 for null terminator
		if (!kmers[i]) {
			fprintf(stderr, "Failed to allocate memory for kmer\n");
			// Free previously allocated memory
			for (int j = 0; j < i; j++) {
				free(kmers[j]);
			}
			free(kmers);
			*size = 0;
			return NULL;
		}
		strncpy(kmers[i], seq + i, k);
		kmers[i][k] = '\0';  // Null terminate
		*size += 1;
	}

	return kmers;
}


// This method counts transcripts and estimates total k-mers
// Returns the number of transcripts
int count_transcripts_and_kmers(char *file, int *total_kmers){
	gzFile fp = gzopen(file, "r");
	if (!fp) {
		fprintf(stderr, "Error: Could not open file '%s'\n", file);
		return -1;
	}
	
	kseq_t *seq = kseq_init(fp);
	int n_transcripts = 0;
	int n_kmers = 0;
	int64_t l;
	
	while ((l = kseq_read(seq)) >= 0) {
		n_transcripts++;
		int seq_len = seq->seq.l;
		if (seq_len >= k) {
			n_kmers += seq_len - k + 1;
		}
	}
	
	kseq_destroy(seq);
	gzclose(fp);
	
	*total_kmers = n_kmers;
	printf("Number of transcripts = %d\n", n_transcripts);
	printf("Estimated total k-mers = %d\n", n_kmers);
	
	return n_transcripts;
}


//store transcript names in t_names to translate final output. The idex is used as the name id.
void store_name_at_idx(char * name, int t_name_id){
	t_names[t_name_id] = calloc(tr_name_chars, sizeof(char));
	strcpy(t_names[t_name_id], name);
}

// stores the k-mer seqence and the label of the transcript it is read from at index kmer_hash(kmer_seq) in the hash table, if the sequence or T-label are not already present
void create_htable(char * tscripts_file)
//void create_htable()
{
	// First pass: count transcripts and k-mers to determine hash table size
	int total_kmers = 0;
	num_t = count_transcripts_and_kmers(tscripts_file, &total_kmers);
	
	if (num_t < 0) {
		fprintf(stderr, "Error: Failed to count transcripts\n");
		exit(1);
	}
	
	// Set hash table size dynamically based on k-mer count
	// Use a load factor of ~0.7 for good performance
	max_hasht_rows = (int)(total_kmers / 0.7);
	if (max_hasht_rows < 1000) max_hasht_rows = 1000;  // Minimum size
	
	printf("Setting hash table size to %d rows (load factor ~0.7)\n", max_hasht_rows);
	
	t_names = calloc(num_t, sizeof(char*));
	t_lengths = calloc(num_t, sizeof(int));
	
	if (!t_names || !t_lengths) {
		fprintf(stderr, "Error: Failed to allocate memory for transcript names/lengths\n");
		exit(1);
	}

	printf("Initializing hash table \n");

	hashTable = calloc(max_hasht_rows, sizeof(kmer_labeled*));
	if (!hashTable) {
		fprintf(stderr, "Error: Failed to allocate memory for hash table\n");
		exit(1);
	}
	
	for(int i = 0; i < max_hasht_rows; i++){
		hashTable[i] = calloc(max_per_row, sizeof(kmer_labeled));
		if (!hashTable[i]) {
			fprintf(stderr, "Error: Failed to allocate memory for hash table row %d\n", i);
			exit(1);
		}
		for(int j = 0; j < max_per_row; j++){
			hashTable[i][j].kmer_seq = calloc(k + 1, sizeof(char));  // +1 for null terminator
			hashTable[i][j].from_reads = 0;
			hashTable[i][j].t_labels = calloc(num_t, sizeof(int)); //every transcript could contain this seq
			if (!hashTable[i][j].kmer_seq || !hashTable[i][j].t_labels) {
				fprintf(stderr, "Error: Failed to allocate memory for hash table entry\n");
				exit(1);
			}
			for(int m = 0; m < num_t; m++){
				hashTable[i][j].t_labels[m] = -1; // -1 means null, there is no label
			}
		}
	}
	printf("Hash table initialized\n");

	// Second pass: read sequences and populate hash table using kseq
	gzFile fp = gzopen(tscripts_file, "r");
	if (!fp) {
		fprintf(stderr, "Error: Could not open transcripts file '%s'\n", tscripts_file);
		exit(1);
	}
	
	kseq_t *seq = kseq_init(fp);
	int t_name_id = 0;
	int64_t l;
	
	while ((l = kseq_read(seq)) >= 0) {
		// Store transcript name
		store_name_at_idx(seq->name.s, t_name_id);
		
		// Store transcript length
		t_lengths[t_name_id] = seq->seq.l;

		// Generate and store k-mers
		int size = 0;
		char** kmers = get_kmers(seq->seq.s, k, &size);
		
		if (!kmers) {
			t_name_id++;
			continue;
		}

		for (int i = 0; i < size; i++){
			char * kmer = kmers[i];
			uint32_t h_row = kmer_hash(kmer);
			int h_col = 0;

			// find matching sequence or find empty seqence (empty if kmer_seq[0] is null from calloc)
			while (strcmp(hashTable[h_row][h_col].kmer_seq, kmer) != 0 && (hashTable[h_row][h_col].kmer_seq[0]) ){
				h_col++;

				if (h_col >= max_per_row){
					fprintf(stderr, "Error: Not enough columns in hash table at row %d: increase max_per_row\n", h_row);
					exit(1);
				}
			}

			strcpy(hashTable[h_row][h_col].kmer_seq, kmer); // extra step if found same seq, rewrites seq

			// add transcript id label at first empty spot in array
			int lab_idx = 0;
			while ((hashTable[h_row][h_col].t_labels[lab_idx] != t_name_id) && (hashTable[h_row][h_col].t_labels[lab_idx] != -1 )){
				lab_idx++;
				if (lab_idx >= num_t){ // should never happen
					fprintf(stderr, "Error: More than max number of transcripts have sequence in common\n");
					lab_idx--;
					break;
				}
			}
			hashTable[h_row][h_col].t_labels[lab_idx] = t_name_id; // extra step if already found, rewrites same int
		}
		
		// Free kmers array
		for (int i = 0; i < size; i++) {
			free(kmers[i]);
		}
		free(kmers);
		t_name_id++;
	}
	
	kseq_destroy(seq);
	gzclose(fp);
	
	printf("Done reading transcripts and storing kmers in hash table!\n\n");
}

void store_read_counts(char * reads_file, int max_reads)
{
	gzFile fp = gzopen(reads_file, "r");
	if (!fp) {
		fprintf(stderr, "Error: Could not open reads file '%s'\n", reads_file);
		exit(1);
	}
	
	kseq_t *seq = kseq_init(fp);
	int64_t l;

	printf("Reading reads \n");

	int reads_processed = 0;
	while ((l = kseq_read(seq)) >= 0 && reads_processed < max_reads) {
		reads_processed++;
		
		if ((int)seq->seq.l < k) {
			continue;  // Skip reads shorter than k-mer size
		}
		
		if (use_full_debruijn) {
			// Full de Bruijn graph approach: check all k-mers in the read
			int seq_len = seq->seq.l;
			int num_kmers = seq_len - k + 1;
			
			// Collect all equivalence classes from matching k-mers
			int **eq_classes = malloc(num_kmers * sizeof(int*));
			if (!eq_classes) {
				continue;
			}
			
			int num_matched = 0;
			char* kmer_buf = malloc((k + 1) * sizeof(char));
			if (!kmer_buf) {
				free(eq_classes);
				continue;
			}
			
			// Try each k-mer in the read
			for (int pos = 0; pos < num_kmers; pos++) {
				strncpy(kmer_buf, seq->seq.s + pos, k);
				kmer_buf[k] = '\0';
				
				// Look up this k-mer in the hash table
				uint32_t h_row = kmer_hash(kmer_buf);
				int h_col = 0;
				
				while (h_col < max_per_row) {
					if (!hashTable[h_row][h_col].kmer_seq[0]) {
						break;
					}
					if (strcmp(hashTable[h_row][h_col].kmer_seq, kmer_buf) == 0) {
						// Found k-mer, store its equivalence class
						int *eq_class = malloc(num_t * sizeof(int));
						if (eq_class) {
							for (int t = 0; t < num_t; t++) {
								eq_class[t] = hashTable[h_row][h_col].t_labels[t];
							}
							eq_classes[num_matched++] = eq_class;
						}
						break;
					}
					h_col++;
				}
			}
			
			free(kmer_buf);
			
			// Compute intersection of all equivalence classes
			if (num_matched > 0) {
				int *intersection = compute_eq_class_intersection(eq_classes, num_matched);
				if (intersection) {
					increment_eq_class_count(intersection);
					free(intersection);
				}
				
				// Free allocated equivalence classes
				for (int i = 0; i < num_matched; i++) {
					free(eq_classes[i]);
				}
			}
			
			free(eq_classes);
			
		} else {
			// Original approach: use first matching k-mer only
			char* read_1kmer = malloc((k + 1) * sizeof(char));
			if (!read_1kmer) {
				fprintf(stderr, "Error: Failed to allocate memory for kmer\n");
				continue;
			}
			strncpy(read_1kmer, seq->seq.s, k);
			read_1kmer[k] = '\0';

			// find read_1kmer in hashTable and increment the from_reads count of it
			uint32_t h_row = kmer_hash(read_1kmer);

			int h_col = 0;
			int k_start = 1;
			int seq_len = seq->seq.l;
			while (strcmp(hashTable[h_row][h_col].kmer_seq, read_1kmer) != 0) {
				if (!(hashTable[h_row][h_col].kmer_seq[0])){
					// Try next kmer in the sequence
					if (k_start + k > seq_len) {
						// No more kmers to try
						break;
					}
					strncpy(read_1kmer, seq->seq.s + k_start, k);
					read_1kmer[k] = '\0';
					k_start++;
					h_row = kmer_hash(read_1kmer);
					h_col = -1;
				}
				h_col++;
			}

			if (hashTable[h_row][h_col].kmer_seq[0]){
				hashTable[h_row][h_col].from_reads++; // only increment if sequence was found in table
			}
			
			free(read_1kmer);
		}
	}
	
	kseq_destroy(seq);
	gzclose(fp);
	
	printf("Done storing read counts (processed %d reads)! \n\n", reads_processed);
}

// Helper function to find matching k-mer in hash table and get equivalence class
// If use_full_debruijn is enabled, computes intersection of all k-mer equivalence classes
// Otherwise, returns the first matching k-mer's equivalence class
int* find_kmer_eq_class(char* seq_str, int seq_len) {
	if (seq_len < k) {
		return NULL;
	}
	
	if (use_full_debruijn) {
		// Full de Bruijn graph: collect all k-mers and compute intersection
		int num_kmers = seq_len - k + 1;
		int **eq_classes = malloc(num_kmers * sizeof(int*));
		if (!eq_classes) {
			return NULL;
		}
		
		int num_matched = 0;
		char* kmer_buf = malloc((k + 1) * sizeof(char));
		if (!kmer_buf) {
			free(eq_classes);
			return NULL;
		}
		
		// Try each k-mer in the read
		for (int pos = 0; pos < num_kmers; pos++) {
			strncpy(kmer_buf, seq_str + pos, k);
			kmer_buf[k] = '\0';
			
			// Look up this k-mer in the hash table
			uint32_t h_row = kmer_hash(kmer_buf);
			int h_col = 0;
			
			while (h_col < max_per_row) {
				if (!hashTable[h_row][h_col].kmer_seq[0]) {
					break;
				}
				if (strcmp(hashTable[h_row][h_col].kmer_seq, kmer_buf) == 0) {
					// Found k-mer, store its equivalence class
					int *eq_class = malloc(num_t * sizeof(int));
					if (eq_class) {
						for (int t = 0; t < num_t; t++) {
							eq_class[t] = hashTable[h_row][h_col].t_labels[t];
						}
						eq_classes[num_matched++] = eq_class;
					}
					break;
				}
				h_col++;
			}
		}
		
		free(kmer_buf);
		
		// Compute intersection of all equivalence classes
		int *intersection = NULL;
		if (num_matched > 0) {
			intersection = compute_eq_class_intersection(eq_classes, num_matched);
			
			// Free allocated equivalence classes
			for (int i = 0; i < num_matched; i++) {
				free(eq_classes[i]);
			}
		}
		
		free(eq_classes);
		return intersection;
		
	} else {
		// Original approach: return first matching k-mer's equivalence class
		char* read_1kmer = malloc((k + 1) * sizeof(char));
		if (!read_1kmer) {
			return NULL;
		}
		
		int k_start = 0;
		int* eq_class = NULL;
		
		// Try multiple k-mers in the sequence
		while (k_start + k <= seq_len) {
			strncpy(read_1kmer, seq_str + k_start, k);
			read_1kmer[k] = '\0';
			
			uint32_t h_row = kmer_hash(read_1kmer);
			int h_col = 0;
			
			while (h_col < max_per_row) {
				if (!hashTable[h_row][h_col].kmer_seq[0]) {
					// Empty slot, k-mer not found in this position
					break;
				}
				if (strcmp(hashTable[h_row][h_col].kmer_seq, read_1kmer) == 0) {
					// Found the k-mer, allocate and copy its equivalence class
					eq_class = malloc(num_t * sizeof(int));
					if (eq_class) {
						for (int t = 0; t < num_t; t++) {
							eq_class[t] = hashTable[h_row][h_col].t_labels[t];
						}
					}
					free(read_1kmer);
					return eq_class;
				}
				h_col++;
			}
			k_start++;
		}
		
		free(read_1kmer);
		return NULL;
	}
}

// Helper function to increment count for first matching k-mer in sequence
int increment_kmer_count(char* seq_str, int seq_len) {
	if (seq_len < k) {
		return 0;
	}
	
	char* read_1kmer = malloc((k + 1) * sizeof(char));
	if (!read_1kmer) {
		return 0;
	}
	
	int k_start = 0;
	int found = 0;
	
	// Try multiple k-mers in the sequence
	while (k_start + k <= seq_len) {
		strncpy(read_1kmer, seq_str + k_start, k);
		read_1kmer[k] = '\0';
		
		uint32_t h_row = kmer_hash(read_1kmer);
		int h_col = 0;
		
		while (h_col < max_per_row) {
			if (!hashTable[h_row][h_col].kmer_seq[0]) {
				// Empty slot, k-mer not found in this position
				break;
			}
			if (strcmp(hashTable[h_row][h_col].kmer_seq, read_1kmer) == 0) {
				// Found the k-mer, increment count
				hashTable[h_row][h_col].from_reads++;
				found = 1;
				free(read_1kmer);
				return 1;
			}
			h_col++;
		}
		k_start++;
	}
	
	free(read_1kmer);
	return found;
}

// Store paired-end read counts
void store_paired_read_counts(char * reads_file1, char * reads_file2, int max_reads)
{
	gzFile fp1 = gzopen(reads_file1, "r");
	gzFile fp2 = gzopen(reads_file2, "r");
	
	if (!fp1) {
		fprintf(stderr, "Error: Could not open reads file '%s'\n", reads_file1);
		exit(1);
	}
	if (!fp2) {
		fprintf(stderr, "Error: Could not open reads file '%s'\n", reads_file2);
		gzclose(fp1);
		exit(1);
	}
	
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);
	int64_t l1, l2;

	printf("Reading paired-end reads \n");

	int reads_processed = 0;
	while (reads_processed < max_reads) {
		l1 = kseq_read(seq1);
		l2 = kseq_read(seq2);
		
		// Check if both reads are available
		if (l1 < 0 || l2 < 0) {
			if (l1 >= 0 || l2 >= 0) {
				fprintf(stderr, "Warning: Unequal number of reads in paired files\n");
			}
			break;
		}
		
		reads_processed++;
		
		// Get equivalence classes for both reads
		int* eq_class1 = find_kmer_eq_class(seq1->seq.s, seq1->seq.l);
		int* eq_class2 = find_kmer_eq_class(seq2->seq.s, seq2->seq.l);
		
		// Combine equivalence classes (intersection of transcripts)
		// If both reads map, use intersection; otherwise use whichever maps
		if (eq_class1 && eq_class2) {
			// Find intersection of both equivalence classes
			int intersection_found = 0;
			for (int i = 0; i < num_t && eq_class1[i] != -1; i++) {
				for (int j = 0; j < num_t && eq_class2[j] != -1; j++) {
					if (eq_class1[i] == eq_class2[j]) {
						intersection_found = 1;
						break;
					}
				}
				if (intersection_found) break;
			}
			
			// If intersection exists, increment counts for the first read's k-mer
			if (intersection_found) {
				increment_kmer_count(seq1->seq.s, seq1->seq.l);
			}
			
			// Free allocated memory
			free(eq_class1);
			free(eq_class2);
		} else if (eq_class1) {
			// Only first read maps
			increment_kmer_count(seq1->seq.s, seq1->seq.l);
			free(eq_class1);
		} else if (eq_class2) {
			// Only second read maps
			increment_kmer_count(seq2->seq.s, seq2->seq.l);
			free(eq_class2);
		}
		// If neither read maps, skip this pair
	}
	
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);
	
	printf("Done storing paired-end read counts (processed %d pairs)! \n\n", reads_processed);
}

// returns the list of T-labels associated with the k-mer sequence from the hash table, including -1s
int * get_equiv_class(char* kmer)
{
	uint32_t h_row = kmer_hash(kmer);
	int h_col = 0;
	while (strcmp(hashTable[h_row][h_col].kmer_seq, kmer) != 0) {
		h_col++;
		if (!(hashTable[h_row][h_col].kmer_seq[0])){
			fprintf(stderr, "kmer sequence %s is not stored in hash table.\n", kmer);
			return NULL;
		}
	}
	// allocate memory to copy transcript ids in eqivalency class from hash table
	int * eq_class = malloc(num_t * sizeof(int));  // Fixed: was sizeof(int*), should be num_t * sizeof(int)
	if (!eq_class) {
		fprintf(stderr, "Failed to allocate memory for equivalence class\n");
		return NULL;
	}
	for (int lab_idx = 0; lab_idx < num_t; lab_idx++){
		eq_class[lab_idx] = hashTable[h_row][h_col].t_labels[lab_idx];
	}
	return eq_class;
}

// Compute intersection of multiple equivalence classes
// Returns a new equivalence class array containing only transcripts present in ALL input classes
int * compute_eq_class_intersection(int **eq_classes, int num_classes) {
	if (num_classes == 0) {
		return NULL;
	}
	
	if (num_classes == 1) {
		// Single class, just copy it
		int *result = malloc(num_t * sizeof(int));
		if (!result) {
			return NULL;
		}
		for (int i = 0; i < num_t; i++) {
			result[i] = eq_classes[0][i];
		}
		return result;
	}
	
	// Start with first class
	int *intersection = malloc(num_t * sizeof(int));
	if (!intersection) {
		return NULL;
	}
	
	// Initialize with all -1
	for (int i = 0; i < num_t; i++) {
		intersection[i] = -1;
	}
	
	int intersection_count = 0;
	
	// For each transcript in the first class
	for (int i = 0; i < num_t && eq_classes[0][i] != -1; i++) {
		int transcript_id = eq_classes[0][i];
		int in_all_classes = 1;
		
		// Check if this transcript is in all other classes
		for (int c = 1; c < num_classes; c++) {
			int found = 0;
			for (int j = 0; j < num_t && eq_classes[c][j] != -1; j++) {
				if (eq_classes[c][j] == transcript_id) {
					found = 1;
					break;
				}
			}
			if (!found) {
				in_all_classes = 0;
				break;
			}
		}
		
		if (in_all_classes) {
			intersection[intersection_count++] = transcript_id;
		}
	}
	
	return intersection;
}

// Increment count for an equivalence class by finding matching k-mer in hash table
void increment_eq_class_count(int *eq_class) {
	if (!eq_class || eq_class[0] == -1) {
		return;
	}
	
	// Find a k-mer in the hash table that has this equivalence class
	// We'll search through the hash table to find a matching entry
	for (int i = 0; i < max_hasht_rows; i++) {
		for (int j = 0; j < max_per_row; j++) {
			if (!hashTable[i][j].kmer_seq[0]) {
				continue;
			}
			
			// Check if this k-mer's equivalence class matches
			int matches = 1;
			for (int t = 0; t < num_t; t++) {
				if (hashTable[i][j].t_labels[t] != eq_class[t]) {
					matches = 0;
					break;
				}
			}
			
			if (matches) {
				hashTable[i][j].from_reads++;
				return;
			}
		}
	}
}


int find_eqc_arr_idx(int * t_labels){

	num_eqc_found++;
	int idx = num_eqc_found; // new class, return if t_labels doesn't match any previously found classes

	for (int i = 0; i < num_eqc_found; i++){
		for (int j = 0; j < num_t; j++){
			if (eqc_arr[i].eq_class_labels[j] != t_labels[j]){ // they are not the same class, try next in table
				break; // breaks out of both loops?
			}
			else if (j == num_t - 1){ // they are the same class since break statement was never reached
				idx = i;
				num_eqc_found--;
			}
		}
	}
	return idx;
}

// extracting the counts of each eq class is used in the EM algorithm, so this will make getting those faster
void store_eqiv_classes()
{
	printf("Initializing equivalency class array \n");
	num_eq_classes = 400;
	printf("num_eq_classes = %d\n", num_eq_classes);

	if (num_eq_classes < 0){
		fprintf(stderr, "Error: Too many transcripts, 2^num_t = %d\n", num_eq_classes);
		exit(1);
	}

	eqc_arr = calloc(num_eq_classes, sizeof(eq_c_count));
	if (!eqc_arr) {
		fprintf(stderr, "Error: Failed to allocate memory for equivalence class array\n");
		exit(1);
	}
	
	for(int i = 0; i < num_eq_classes; i++){
		eqc_arr[i].count = 0;
		eqc_arr[i].eq_class_id = -1;
		eqc_arr[i].eq_class_labels = calloc(num_t, sizeof(int));
		if (!eqc_arr[i].eq_class_labels) {
			fprintf(stderr, "Error: Failed to allocate memory for equivalence class labels\n");
			exit(1);
		}
		for(int j = 0; j < num_t; j++){
			eqc_arr[i].eq_class_labels[j] = -1; // important for comparing equivalency classes
		}
	}

	for(int i = 0; i < max_hasht_rows; i++){
		for(int j = 0; j < max_per_row; j++){
			if (hashTable[i][j].from_reads > 0){
				int idx = -1;
				idx = find_eqc_arr_idx(hashTable[i][j].t_labels); // find the equivalency class of the kmer
				eqc_arr[idx].count += hashTable[i][j].from_reads;

				if (idx == num_eqc_found){
					eqc_arr[idx].eq_class_id = idx; // changes from -1 if new equivalency class
					for(int m = 0; m < num_t; m++){
						eqc_arr[idx].eq_class_labels[m] = hashTable[i][j].t_labels[m];
					}
				}

			}

		}
	}
}

// This method print all the eqivalence classes that were found
void print_eqc_arr(){
	int col_print = 15;

	for (int i = 0 ; i < num_eqc_found; i++){
		printf("eqc_arr[%2d].eq_class_labels = ", i);
		for (int j = 0; j < col_print; j++){
			//if(eqc_arr[i].eq_class_labels[j] != -1)
			printf("%4d ", eqc_arr[i].eq_class_labels[j]);
		}
		printf("count = %6d \n", eqc_arr[i].count);
	}
}

// This method prints the hash table and it's content
void print_hashTable(){
	int col_print = 15;
	for (int i = 0 ; i < max_hasht_rows; i++){
		for (int j = 0; j < max_per_row; j++){
			if (hashTable[i][j].kmer_seq[0]){
				printf("hashTable[%3d][%3d].kmer_seq = %s\n", i, j, hashTable[i][j].kmer_seq);
				printf("hashTable[%3d][%3d].from_reads = %d\n", i, j, hashTable[i][j].from_reads);
				printf("hashTable[%3d][%3d].t_labels = ", i, j); //every transcript could contain this seq
				for(int m = 0; m < col_print; m++){
					printf("%4d ", hashTable[i][j].t_labels[m]); // -1 means null, there is no label
				}
				printf("\n");
			}
		}
	}
}

// Here we are printing the vector of lengths
void print_t_lengths(){
	for (int i = 0; i < num_t; i++){  // Fixed: was <=, should be <
		printf("t_lengths[%d] = %d\n", i, t_lengths[i]);
	}
}


// Cleans the given text by changing the next line symbol for a 
// string end symbol
char *clean(char *str){

	for(int i=0; str[i]!='\0'; i++){
		if(str[i]=='\n')
			str[i]='\0';
	}

	return str;

}


// This is an implementation of fast hash table that was used in the paper for kallisto
// from https://github.com/jwerle/murmurhash.c
uint32_t murmurhash (const char *key, uint32_t len, uint32_t seed) {
	uint32_t c1 = 0xcc9e2d51;
	uint32_t c2 = 0x1b873593;
	uint32_t r1 = 15;
	uint32_t r2 = 13;
	uint32_t m = 5;
	uint32_t n = 0xe6546b64;
	uint32_t h = 0;
	uint32_t k = 0;
	uint8_t *d = (uint8_t *) key; // 32 bit extract from `key'
	const uint32_t *chunks = NULL;
	const uint8_t *tail = NULL; // tail - last 8 bytes
	int i = 0;
	int l = len / 4; // chunk length

	h = seed;

	chunks = (const uint32_t *) (d + l * 4); // body
	tail = (const uint8_t *) (d + l * 4); // last 8 byte chunk of `key'

	// for each 4 byte chunk of `key'
	for (i = -l; i != 0; ++i) {
		// next 4 byte chunk of `key'
		k = chunks[i];

		// encode next 4 byte chunk of `key'
		k *= c1;
		k = (k << r1) | (k >> (32 - r1));
		k *= c2;

		// append to hash
		h ^= k;
		h = (h << r2) | (h >> (32 - r2));
		h = h * m + n;
	}

	k = 0;

	// remainder
	switch (len & 3) { // `len % 4'
		case 3: 
			k ^= (tail[2] << 16);
			/* fall through */
		case 2: 
			k ^= (tail[1] << 8);
			/* fall through */
		case 1:
			k ^= tail[0];
			k *= c1;
			k = (k << r1) | (k >> (32 - r1));
			k *= c2;
			h ^= k;
	}

	h ^= len;

	h ^= (h >> 16);
	h *= 0x85ebca6b;
	h ^= (h >> 13);
	h *= 0xc2b2ae35;
	h ^= (h >> 16);

	return h;
}


// Parallel worker function for update_alphas
static void worker_update_alphas(void *data, long k, int tid) {
	(void)tid; // Unused parameter
	em_data_t *em_data = (em_data_t*)data;
	double sum = 0;
	
	for(int i = 1; i < num_eqc_found; i++){
		// KSD correction: multiply by count
		sum += eqc_arr[i].count * em_data->transcript_prob[i][k];
	}
	
	em_data->alpha[k] = sum/num_eqc_found;
}

// Parallel worker function for update_trans_probs
static void worker_update_trans_probs(void *data, long idx, int tid) {
	(void)tid; // Unused parameter
	em_data_t *em_data = (em_data_t*)data;
	long i = idx + 1; // Offset by 1 since we start from 1
	int j = 0;
	double sum = 0;
	
	// Iterate the transcripts
	// Here we are adding the prbability of a specific equivalence class
	while(eqc_arr[i].eq_class_labels[j] != -1){
		sum += em_data->probs[eqc_arr[i].eq_class_labels[j]];
		j++;
	}
	
	for(int k = 0; k < j; k++)
		em_data->transcript_prob[i][eqc_arr[i].eq_class_labels[k]] = em_data->probs[eqc_arr[i].eq_class_labels[k]] / sum;
}

// This is our implementation of the EM algorithm
void EM(int *lengths, double eps){
	if (bootstrap_samples > 0) {
		bootstrap_EM(lengths, eps);
	} else {
		// Regular EM without bootstrap
		double *probs = calloc(num_t + 1, sizeof(double));
		if (!probs) {
			fprintf(stderr, "Error: Failed to allocate memory for EM algorithm\n");
			exit(1);
		}
		
		run_EM_iteration(lengths, eps, probs);
		
		// Write abundances to output file
		if (g_output_file) {
			write_abundances(g_output_file, probs);
		}
		
		// This is the output of the EM (for validation if real.txt exists)
		compare(probs);
		
		// Free the probs we allocated here
		free(probs);
	}
}

// Run a single EM iteration (used for both regular and bootstrap)
void run_EM_iteration(int *lengths, double eps, double *probs_out) {
	// Because the first cell is null in the eqc_arr
	double *probs = calloc(num_t + 1, sizeof(double));
	double *alpha = calloc(num_t + 1,  sizeof(double));
	
	if (!probs || !alpha) {
		fprintf(stderr, "Error: Failed to allocate memory for EM algorithm\n");
		exit(1);
	}

	double llk = 0;
	double prev_llk = 0;
	double **transcirpt_prob = alloc_matrix(num_eqc_found, num_t + 1);
	
	if (!transcirpt_prob) {
		fprintf(stderr, "Error: Failed to allocate memory for transcript probabilities\n");
		free(probs);
		free(alpha);
		exit(1);
	}

	// Initialize the probs
	for(int i = 0; i <= num_t; i++){
		probs[i] = 1./(num_t + 1);
		alpha[i] = 0;
	}

	double diff;
	int iterations = 0;
	do {
		prev_llk = llk;
		update_parameters(alpha, probs, lengths, transcirpt_prob);
		llk = likelihood(alpha, lengths);
		diff = fabs(llk - prev_llk);
		iterations++;
		if (iterations % 10 == 0 && bootstrap_samples == 0) {
			fprintf(stderr, "Iteration %d: llk = %f, diff = %e\n", iterations, llk, diff);
		}
	} while(diff > eps);
	
	if (bootstrap_samples == 0) {
		printf("\nEM converged after %d iterations\n", iterations);
		printf("Final log-likelihood: %f\n", llk);
	}
	
	// Copy results to output
	for(int i = 0; i <= num_t; i++){
		probs_out[i] = probs[i];
	}
	
	// Free temporary memory
	for(int i = 0; i < num_eqc_found; i++){
		free(transcirpt_prob[i]);
	}
	free(transcirpt_prob);
	free(probs);
	free(alpha);
}

// Bootstrap EM algorithm
void bootstrap_EM(int *lengths, double eps) {
	printf("Running bootstrap with %d samples...\n", bootstrap_samples);
	
	// Allocate memory for bootstrap results
	double **bootstrap_results = alloc_matrix(bootstrap_samples, num_t + 1);
	if (!bootstrap_results) {
		fprintf(stderr, "Error: Failed to allocate memory for bootstrap results\n");
		exit(1);
	}
	
	// Save original counts
	int *original_counts = malloc(num_eqc_found * sizeof(int));
	if (!original_counts) {
		fprintf(stderr, "Error: Failed to allocate memory for original counts\n");
		exit(1);
	}
	for(int i = 0; i < num_eqc_found; i++){
		original_counts[i] = eqc_arr[i].count;
	}
	
	// Run original EM first
	double *original_probs = calloc(num_t + 1, sizeof(double));
	if (!original_probs) {
		fprintf(stderr, "Error: Failed to allocate memory for original probs\n");
		free(original_counts);
		exit(1);
	}
	
	printf("Running original EM...\n");
	run_EM_iteration(lengths, eps, original_probs);
	printf("EM converged\n");
	
	// Run bootstrap iterations
	for(int bs = 0; bs < bootstrap_samples; bs++) {
		if ((bs + 1) % 10 == 0 || bs == 0) {
			printf("Bootstrap iteration %d/%d\n", bs + 1, bootstrap_samples);
		}
		
		// Proper bootstrap resampling with replacement
		// First, calculate total read count
		int total_reads = 0;
		for(int i = 1; i < num_eqc_found; i++){
			total_reads += original_counts[i];
		}
		
		// Reset all counts to 0
		for(int i = 1; i < num_eqc_found; i++){
			eqc_arr[i].count = 0;
		}
		
		// Resample total_reads times with replacement
		for(int r = 0; r < total_reads; r++){
			// Use rejection sampling for uniform distribution
			int random_idx;
			int range = num_eqc_found - 1;
			int limit = RAND_MAX - (RAND_MAX % range);
			do {
				random_idx = rand();
			} while (random_idx >= limit);
			random_idx = 1 + (random_idx % range);
			
			eqc_arr[random_idx].count++;
		}
		
		// Run EM on resampled data
		run_EM_iteration(lengths, eps, bootstrap_results[bs]);
	}
	
	// Restore original counts
	for(int i = 0; i < num_eqc_found; i++){
		eqc_arr[i].count = original_counts[i];
	}
	
	printf("Bootstrap complete. Computing statistics...\n");
	
	// Write results with bootstrap statistics
	if (g_output_file) {
		write_abundances_with_bootstrap(g_output_file, original_probs, bootstrap_results, bootstrap_samples);
	}
	
	// This is the output of the EM (for validation if real.txt exists)
	compare(original_probs);
	
	// Free memory
	free(original_probs);
	free(original_counts);
	for(int i = 0; i < bootstrap_samples; i++){
		free(bootstrap_results[i]);
	}
	free(bootstrap_results);
}

// This methods update the alphas for each of the transcripts
void update_alphas(double *alpha, double **transcirpt_prob){
	if (n_threads > 1) {
		em_data_t em_data = {alpha, transcirpt_prob, NULL, NULL};
		kt_for(n_threads, worker_update_alphas, &em_data, num_t + 1);
	} else {
		for(int k = 0; k <= num_t; k++){
			double sum = 0;

			for(int i = 1; i < num_eqc_found; i++){
				// KSD correction: multiply by count
				sum += eqc_arr[i].count * transcirpt_prob[i][k];
			}

			alpha[k] = sum/num_eqc_found;
		}
	}
}

// This updates the probabilities
void update_prob(double *alpha, double *probs, int *lengths){

	double sum = 0;
	for(int j = 0; j <= num_t; j++){
		sum += alpha[j] / lengths[j];
	}

	for(int k = 0; k <= num_t; k++){
		probs[k] = (alpha[k] / lengths[k]) / sum;
	}

}

// Update the transcriptions probabilities
void update_trans_probs(double *probs, double **transcirpt_prob){
	if (n_threads > 1) {
		em_data_t em_data = {NULL, transcirpt_prob, probs, NULL};
		// Start from 1 because first cell is null
		kt_for(n_threads, worker_update_trans_probs, &em_data, num_eqc_found - 1);
		// Note: kt_for processes indices 0 to n-1, but we want 1 to num_eqc_found-1
		// So we need to adjust the worker function
	} else {
		for(int i = 1; i < num_eqc_found; i++){
			int j = 0;
			double sum = 0;

			// Iterate the transcripts
			// Here we are adding the prbability of a specific equivalence class
			while(eqc_arr[i].eq_class_labels[j] != -1){
				sum += probs[eqc_arr[i].eq_class_labels[j]];
				j++;
			}

			for(int k = 0; k < j; k++)
				transcirpt_prob[i][eqc_arr[i].eq_class_labels[k]] = probs[eqc_arr[i].eq_class_labels[k]] / sum;
		}
	}
}

// Here we are computing the likelihood to determine convergence
double likelihood(double *alpha, int *lengths){
	double llk = 0;

	// Iterate through the equivalence class vector
	for(int e = 1; e < num_eqc_found; e++){
		double sum = 0;
		int t = 0;

		// Iterate through the transcript in a class
		while(eqc_arr[e].eq_class_labels[t] != -1){
			sum += (alpha[eqc_arr[e].eq_class_labels[t]]) / (lengths[eqc_arr[e].eq_class_labels[t]]);
			t++;
		}
		llk += eqc_arr[e].count * log(sum);
	}

	return llk;
}

// From here we call all the update methods
void update_parameters(double *alpha, double *probs, int *lengths, double **transcirpt_prob)
{
	// First, we calculate the probabilities of the trascripts
	update_trans_probs(probs, transcirpt_prob);

	// Second, we update the alphas
	update_alphas(alpha, transcirpt_prob);

	// We update the class probabilities
	update_prob(alpha, probs, lengths);
}

// This method allocate a 2D array of size n and r in the memory
double **alloc_matrix(int r, int c){
	double **mat = calloc(r, sizeof(double*));
	
	if (!mat) {
		fprintf(stderr, "Error: Failed to allocate memory for matrix\n");
		return NULL;
	}

	for(int i = 0; i < r; i++) {
		mat[i] = calloc(c, sizeof(double));
		if (!mat[i]) {
			fprintf(stderr, "Error: Failed to allocate memory for matrix row %d\n", i);
			// Free previously allocated rows
			for (int j = 0; j < i; j++) {
				free(mat[j]);
			}
			free(mat);
			return NULL;
		}
	}

	return mat;
}

// Here we free all the spaced used by the hash table
void free_hashTable(){
	for(int i = 0; i < max_hasht_rows; i++){
		for(int j = 0; j < max_per_row; j++){
			free(hashTable[i][j].kmer_seq);
			free(hashTable[i][j].t_labels);
		}
		free(hashTable[i]);
	}
	free(hashTable);
}

// This method calculate the mean squared error
void compare(double *probs){

	FILE * tscripts_filep = fopen("real.txt", "r");
	if (!tscripts_filep) {
		fprintf(stderr, "Warning: Could not open real.txt for comparison\n");
		return;
	}
	
	char n[15];
	int d1;
	int *d = malloc(num_t * sizeof(int));
	if (!d) {
		fprintf(stderr, "Error: Failed to allocate memory for comparison\n");
		fclose(tscripts_filep);
		return;
	}
	
	int tot = 0;

	for(int i = 0; i < num_t; i++){
		if (fscanf(tscripts_filep, "%s %d %d", n, &d1, &d[i]) != 3) {
			fprintf(stderr, "Warning: Error reading real.txt at line %d\n", i);
			break;
		}
		tot += d[i];
	}

	//Compare the real and the predicted
	double error = 0;
	printf("Labels \t\t Real \t\t Predicted\n");
	for(int i = 0; i < num_t; i++){
		double p = (double) d[i] / tot;
		error += pow((p - probs[i + 1]), 2);
		printf("%s \t %lf \t %lf\n", t_names[i], p, probs[i + 1]);  // Fixed: t_names is 0-indexed
	}

	printf("The MSE is %lf\n", (error/num_t));

	free(d);
	fclose(tscripts_filep);
}

// Write transcript abundances to output file in TSV format
void write_abundances(const char *output_file, double *probs) {
	FILE *out = fopen(output_file, "w");
	if (!out) {
		fprintf(stderr, "Error: Could not open output file '%s'\n", output_file);
		return;
	}
	
	// Write header
	fprintf(out, "target_id\tlength\teff_length\test_counts\ttpm\n");
	
	// Calculate total TPM
	double total_tpm = 0;
	for (int i = 1; i <= num_t; i++) {
		total_tpm += probs[i];
	}
	
	// Write each transcript
	for (int i = 1; i <= num_t; i++) {
		double tpm = (total_tpm > 0) ? (probs[i] / total_tpm) * 1000000.0 : 0;
		double est_counts = probs[i] * t_lengths[i - 1];
		
		fprintf(out, "%s\t%d\t%d\t%.2f\t%.4f\n", 
			t_names[i - 1],  // Fixed: t_names is 0-indexed
			t_lengths[i - 1],
			t_lengths[i - 1],  // effective length (simplified, same as length)
			est_counts,
			tpm);
	}
	
	fclose(out);
	printf("Abundances written to %s\n", output_file);
}

// Write transcript abundances with bootstrap statistics to output file in TSV format
void write_abundances_with_bootstrap(const char *output_file, double *probs, double **bootstrap_results, int n_bootstrap) {
	FILE *out = fopen(output_file, "w");
	if (!out) {
		fprintf(stderr, "Error: Could not open output file '%s'\n", output_file);
		return;
	}
	
	// Write header with bootstrap columns
	fprintf(out, "target_id\tlength\teff_length\test_counts\ttpm\tbs_mean_tpm\tbs_std_tpm\n");
	
	// Calculate total TPM for original
	double total_tpm = 0;
	for (int i = 1; i <= num_t; i++) {
		total_tpm += probs[i];
	}
	
	// Write each transcript with bootstrap statistics
	for (int i = 1; i <= num_t; i++) {
		double tpm = (total_tpm > 0) ? (probs[i] / total_tpm) * 1000000.0 : 0;
		double est_counts = probs[i] * t_lengths[i - 1];
		
		// Calculate bootstrap mean and std dev for this transcript
		double bs_sum = 0;
		double bs_sum_sq = 0;
		for (int bs = 0; bs < n_bootstrap; bs++) {
			// Calculate TPM for this bootstrap sample
			double bs_total_tpm = 0;
			for (int j = 1; j <= num_t; j++) {
				bs_total_tpm += bootstrap_results[bs][j];
			}
			double bs_tpm = (bs_total_tpm > 0) ? (bootstrap_results[bs][i] / bs_total_tpm) * 1000000.0 : 0;
			bs_sum += bs_tpm;
			bs_sum_sq += bs_tpm * bs_tpm;
		}
		double bs_mean = bs_sum / n_bootstrap;
		double bs_variance = (bs_sum_sq / n_bootstrap) - (bs_mean * bs_mean);
		double bs_std = sqrt(bs_variance > 0 ? bs_variance : 0);
		
		fprintf(out, "%s\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\n", 
			t_names[i - 1],  // t_names is 0-indexed
			t_lengths[i - 1],
			t_lengths[i - 1],  // effective length (simplified, same as length)
			est_counts,
			tpm,
			bs_mean,
			bs_std);
	}
	
	fclose(out);
	printf("Abundances with bootstrap statistics written to %s\n", output_file);
}


