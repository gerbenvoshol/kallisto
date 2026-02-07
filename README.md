# Kallisto C Implementation

A reimplementation of kallisto, an RNA-seq quantification tool, in C. This is a lightweight version that implements the core quantification algorithm using k-mer indexing and Expectation-Maximization (EM) for transcript abundance estimation.

## Features

- **K-mer Indexing**: Builds a hash table index from transcriptome FASTA/FASTQ files
- **FASTQ/FASTA Support**: Uses kseq.h for efficient reading of both FASTA and FASTQ formats
- **Gzipped File Support**: Can read compressed .gz files directly
- **Dynamic Hash Table**: Hash table size automatically scales based on k-mer count
- **Pseudoalignment**: Maps reads to transcripts using k-mer matching
- **Paired-end Read Support**: Process paired-end reads for improved accuracy
- **Equivalence Classes**: Groups reads by their compatible transcripts
- **EM Algorithm**: Estimates transcript abundances using the EM algorithm
- **Bootstrap Analysis**: Statistical uncertainty estimation through bootstrapping
- **Multi-threading Support**: Uses kthread.h for parallel computation in EM algorithm
- **Command-line Interface**: Full argument parsing with help and version info
- **TSV Output**: Produces standard kallisto-compatible output format

## Building

### Prerequisites
- GCC compiler (or compatible C compiler)
- zlib development library (for gzipped file support)
- pthread library (for multi-threading support)
- Make (optional, for using Makefile)

On Ubuntu/Debian:
```bash
sudo apt-get install build-essential zlib1g-dev
```

On macOS:
```bash
brew install gcc zlib
```

### Compilation

Using Make:
```bash
make
```

Or manually:
```bash
gcc -Wall -Wextra -O2 -o kallisto kallisto.c kthread.c -lm -lz -lpthread
gcc -Wall -Wextra -O2 -o binary_convert Binary_convert.c
```

## Usage

### Basic Usage

```bash
# Single-end mode with FASTA format
./kallisto -i transcriptome.fasta -r reads.fasta -o output.tsv

# Single-end mode with FASTQ format
./kallisto -i transcriptome.fastq -r reads.fastq -o output.tsv

# Paired-end mode
./kallisto -i transcriptome.fasta -1 reads_1.fastq -2 reads_2.fastq -o output.tsv

# Gzipped files (automatically detected)
./kallisto -i transcriptome.fasta.gz -r reads.fastq.gz -o output.tsv

# Paired-end with gzipped files
./kallisto -i transcriptome.fasta.gz -1 reads_1.fastq.gz -2 reads_2.fastq.gz -o output.tsv

# Use multiple threads for faster processing
./kallisto -i transcriptome.fasta -r reads.fasta -o output.tsv -t 4

# With bootstrap for uncertainty estimation
./kallisto -i transcriptome.fasta -r reads.fasta -o output.tsv -b 100
```

### Command-line Options

```
Required arguments:
  -i, --index <file>        Transcriptome index/reference file (FASTA/FASTQ format)
  -o, --output <file>       Output file for abundance estimates

Single-end mode:
  -r, --reads <file>        Input reads file (FASTA/FASTQ format, gzipped supported)

Paired-end mode:
  -1, --reads1 <file>       First reads file (FASTA/FASTQ format, gzipped supported)
  -2, --reads2 <file>       Second reads file (FASTA/FASTQ format, gzipped supported)

Optional arguments:
  -k, --kmer-size <int>     K-mer size (default: 31)
  -e, --epsilon <float>     EM convergence threshold (default: 0.01)
  -m, --max-reads <int>     Maximum number of reads to process (default: all)
  -t, --threads <int>       Number of threads (default: 1)
  -b, --bootstrap <int>     Number of bootstrap samples (default: 0)
  -d, --full-debruijn       Use full de Bruijn graph (all k-mers) for pseudoalignment
  -h, --help                Display help message
  -v, --version             Display version information
```

### Example

```bash
# Quantify transcript abundances with FASTA (single-end)
./kallisto -i transcriptome.fasta -r sample.fasta -o abundances.tsv -k 31 -e 0.01

# Quantify with FASTQ format (single-end)
./kallisto -i transcriptome.fasta -r sample.fastq -o abundances.tsv -k 31 -e 0.01

# Quantify with paired-end reads
./kallisto -i transcriptome.fasta -1 sample_1.fastq -2 sample_2.fastq -o abundances.tsv -k 31 -e 0.01

# Quantify with gzipped files
./kallisto -i transcriptome.fasta.gz -r sample.fastq.gz -o abundances.tsv -k 31 -e 0.01

# Quantify with bootstrap (100 samples)
./kallisto -i transcriptome.fasta -r sample.fastq -o abundances.tsv -b 100

# Paired-end with bootstrap and multiple threads
./kallisto -i transcriptome.fasta -1 sample_1.fastq.gz -2 sample_2.fastq.gz -o abundances.tsv -b 100 -t 4

# With full de Bruijn graph for more accurate pseudoalignment
./kallisto -i transcriptome.fasta -r sample.fastq -o abundances.tsv -d

# Paired-end with full de Bruijn graph
./kallisto -i transcriptome.fasta -1 sample_1.fastq.gz -2 sample_2.fastq.gz -o abundances.tsv -d

# Convert FASTA file (if needed)
./binary_convert input.fasta output.bin
```

## Input Format

### Transcriptome File
Both FASTA and FASTQ formats are supported. Files can be gzipped (.gz extension).

**FASTA format:**
```
>transcript1
ATCGATCGATCG...
>transcript2
GCTAGCTAGCTA...
```

**FASTQ format:**
```
@transcript1
ATCGATCGATCG...
+
IIIIIIIIIIII...
@transcript2
GCTAGCTAGCTA...
+
IIIIIIIIIIII...
```

### Reads File
Both FASTA and FASTQ formats are supported. Files can be gzipped (.gz extension).

**FASTA format:**
```
>read1
ATCGATCGATCG...
>read2
GCTAGCTAGCTA...
```

**FASTQ format:**
```
@read1
ATCGATCGATCG...
+
IIIIIIIIIIII...
@read2
GCTAGCTAGCTA...
+
IIIIIIIIIIII...
```

## Output Format

The output is a tab-separated values (TSV) file with the following columns:

### Standard Output (without bootstrap)
- `target_id`: Transcript identifier
- `length`: Transcript length
- `eff_length`: Effective length (currently same as length)
- `est_counts`: Estimated counts
- `tpm`: Transcripts per million

Example output:
```
target_id    length    eff_length    est_counts    tpm
transcript1  1000      1000          150.00        75000.0000
transcript2  1500      1500          200.00        66666.6667
```

### Bootstrap Output (with -b flag)
When bootstrap is enabled, additional columns are added:
- `target_id`: Transcript identifier
- `length`: Transcript length
- `eff_length`: Effective length (currently same as length)
- `est_counts`: Estimated counts
- `tpm`: Transcripts per million
- `bs_mean_tpm`: Bootstrap mean TPM
- `bs_std_tpm`: Bootstrap standard deviation of TPM

Example output with bootstrap:
```
target_id    length    eff_length    est_counts    tpm         bs_mean_tpm  bs_std_tpm
transcript1  1000      1000          150.00        75000.0000  74500.0000   1200.5000
transcript2  1500      1500          200.00        66666.6667  66800.0000   800.3000
```

## Algorithm Overview

1. **Indexing**: Build k-mer hash table from transcriptome
2. **Pseudoalignment**: Map reads to transcripts using k-mers
3. **Equivalence Classes**: Group reads by compatible transcript sets
4. **EM Algorithm**: Iteratively estimate transcript abundances until convergence

## Implementation Notes

### Improvements from Original
- **Integrated kseq.h**: Efficient reading of FASTA/FASTQ formats
- **FASTQ support**: Can now process FASTQ files in addition to FASTA
- **Gzipped file support**: Direct reading of .gz compressed files
- **Dynamic hash table sizing**: Hash table size automatically scales based on k-mer count (load factor ~0.7)
- **Paired-end read support**: Process paired-end reads for improved accuracy
- **Bootstrap analysis**: Statistical uncertainty estimation with configurable sample count
- **Multi-threading support**: Integrated kthread.h for parallel computation in EM algorithm
  - Parallelized update_alphas function across transcripts
  - Parallelized update_trans_probs function across equivalence classes
  - Thread count configurable via -t/--threads option
- Fixed critical memory bugs (boundary conditions, allocation errors)
- Added comprehensive error checking for all file operations
- Implemented proper memory management (fixed memory leaks)
- Added command-line argument parsing
- Added standard TSV output format
- Improved k-mer size calculation (off-by-one fix)
- Added null termination for all strings
- Fixed compiler warnings

### Limitations
- Simplified effective length calculation
- Paired-end support uses simplified intersection approach

### Pseudoalignment Options
The implementation supports two pseudoalignment modes:
1. **First k-mer only (default)**: Uses only the first matching k-mer from each read for pseudoalignment. This is faster but less accurate.
2. **Full de Bruijn graph (--full-debruijn flag)**: Examines all k-mers in each read and computes the intersection of their equivalence classes. This is more accurate and matches the approach used by the original kallisto, but is computationally more intensive.

### Performance Considerations
- Hash table size dynamically scales with k-mer count (load factor ~0.7)
- Memory usage: O(k * n) where k = k-mer size, n = number of transcripts
- Time complexity: O(m * k) for indexing + O(r * k) for quantification
  where m = total transcript length, r = number of reads
- kseq.h provides efficient buffered I/O for reading FASTA/FASTQ files

## Technical Details

### Hash Function
Uses MurmurHash3 for k-mer hashing, same as the original kallisto implementation.

### EM Convergence
Convergence is determined when the change in log-likelihood between iterations falls below epsilon (default: 0.01).

### Data Structures
- **K-mer Hash Table**: 2D array with collision handling
- **Equivalence Classes**: Array of counts and transcript label sets
- **Transcript Probability Matrix**: 2D array for EM computation

## Development

### Building in Debug Mode
```bash
make debug
```

### Running Tests
Note: Test data is not included in this repository. You can use your own transcriptome and reads files.

## Project History

Originally developed as a group project for STAT590. This version includes numerous bug fixes, improvements, and feature additions to approach feature parity with the original C++ kallisto implementation.

## References

- Original kallisto paper: Bray et al., Nature Biotechnology (2016)
- kallisto GitHub: https://github.com/pachterlab/kallisto

## License

This is an educational reimplementation. Please refer to the original kallisto license for production use.

## Authors

- Original implementation: STAT590 group project
- Bug fixes and improvements: Multiple contributors
- Current maintainer: See git history

## Acknowledgments

- Original kallisto authors (Páll Melsted, Harold Pimentel, Lior Pachter)
- MurmurHash3 implementation by Austin Appleby
