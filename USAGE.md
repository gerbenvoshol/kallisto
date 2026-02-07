# Kallisto C Implementation - Usage Examples

## Quick Start

### 1. Build the Programs
```bash
make
```

This creates two executables:
- `kallisto` - The main quantification program
- `binary_convert` - FASTA file converter utility

### 2. Prepare Your Data

You need two FASTA files:
- **Transcriptome reference**: Contains all transcript sequences
- **Reads file**: Contains RNA-seq read sequences

### 3. Run Quantification

```bash
./kallisto -i transcriptome.fasta -r reads.fasta -o abundances.tsv
```

## Detailed Examples

### Example 1: Basic Quantification

```bash
# With default k-mer size (31)
./kallisto -i reference.fasta -r sample1.fasta -o sample1_abundance.tsv

# With custom k-mer size
./kallisto -i reference.fasta -r sample1.fasta -o sample1_abundance.tsv -k 25

# With stricter convergence
./kallisto -i reference.fasta -r sample1.fasta -o sample1_abundance.tsv -e 0.001
```

### Example 2: Processing Multiple Samples

```bash
# Process multiple samples in a loop
for sample in sample1 sample2 sample3; do
    ./kallisto -i reference.fasta \
               -r ${sample}.fasta \
               -o ${sample}_abundance.tsv \
               -k 31 \
               -e 0.01
done
```

### Example 3: Limited Read Processing

```bash
# Process only first 10000 reads (for testing)
./kallisto -i reference.fasta -r reads.fasta -o test_output.tsv -m 10000
```

### Example 4: Using Binary Convert

```bash
# Convert FASTA to binary format (if needed)
./binary_convert input.fasta output.bin

# Then use the binary file
./kallisto -i output.bin -r reads.fasta -o abundances.tsv
```

## Understanding the Output

The output TSV file contains the following columns:

1. **target_id**: Transcript identifier (from FASTA header)
2. **length**: Transcript length in nucleotides
3. **eff_length**: Effective length (currently same as length)
4. **est_counts**: Estimated read counts for this transcript
5. **tpm**: Transcripts Per Million (normalized abundance)

Example output:
```
target_id       length  eff_length  est_counts  tpm
transcript_A    1500    1500        234.56      156789.1234
transcript_B    2000    2000        123.45      61725.0000
```

## Common Use Cases

### Differential Expression Analysis Setup

```bash
# Quantify control samples
for i in {1..3}; do
    ./kallisto -i reference.fasta \
               -r control_${i}.fasta \
               -o control_${i}_abundance.tsv
done

# Quantify treatment samples
for i in {1..3}; do
    ./kallisto -i reference.fasta \
               -r treatment_${i}.fasta \
               -o treatment_${i}_abundance.tsv
done
```

### Quality Control Run

```bash
# Test with small subset
./kallisto -i reference.fasta \
           -r reads.fasta \
           -o qc_test.tsv \
           -m 1000 \
           -k 25
```

## Parameter Selection Guide

### K-mer Size (`-k`)
- **Default: 31** - Good for most applications
- **Smaller (21-25)**: Better for shorter reads or lower quality data
- **Larger (35-51)**: Better for longer reads and high quality data
- Must be odd number for biological reasons
- Should be shorter than read length

### Epsilon (`-e`)
- **Default: 0.01** - Good balance of speed and accuracy
- **Smaller (0.001-0.005)**: More accurate but slower
- **Larger (0.05-0.1)**: Faster but less accurate

### Max Reads (`-m`)
- **Default: unlimited** - Process all reads
- Use smaller values for:
  - Testing/debugging
  - Quick quality checks
  - Memory-constrained systems

## Troubleshooting

### Out of Memory
- Reduce k-mer size
- Process fewer reads with `-m`
- Increase system memory

### Slow Convergence
- Increase epsilon slightly
- Check for low-quality data
- Verify transcriptome reference quality

### No Reads Mapping
- Check k-mer size compatibility
- Verify FASTA format correctness
- Ensure transcriptome matches organism

### NaN in Output
- May occur with very small test datasets
- Verify sufficient read depth
- Check for numerical stability issues

## Best Practices

1. **Choose appropriate k-mer size**
   - For 150bp reads: k=31 (default)
   - For 100bp reads: k=25
   - For 50bp reads: k=21

2. **Use consistent parameters**
   - Keep k-mer size constant across samples in an experiment
   - Use same reference for all comparisons

3. **Validate results**
   - Check that TPM values sum to ~1,000,000
   - Verify expected highly-expressed genes appear
   - Compare with biological expectations

4. **Performance optimization**
   - Use binary format for repeated analyses
   - Process samples in parallel on multi-core systems
   - Consider limiting reads during development/testing

## Integration with Downstream Analysis

### Import into R
```r
# Read abundance file
abundance <- read.table("abundances.tsv", header=TRUE, sep="\t")

# Calculate total counts
total_counts <- sum(abundance$est_counts)

# Filter low abundance
filtered <- abundance[abundance$tpm > 1, ]
```

### Import into Python
```python
import pandas as pd

# Read abundance file
df = pd.read_csv('abundances.tsv', sep='\t')

# Filter and analyze
high_expressed = df[df['tpm'] > 100]
print(f"Found {len(high_expressed)} highly expressed transcripts")
```

## Additional Resources

- See README.md for algorithm details
- Check kallisto_project.pdf for theoretical background
- Original kallisto: https://pachterlab.github.io/kallisto/

## Getting Help

```bash
./kallisto --help     # Display help message
./kallisto --version  # Show version information
```

For bugs or issues, please report with:
- Command used
- Sample data characteristics (number of transcripts, reads)
- Error messages or unexpected output
