# FunSeq2 Python Conversion

This directory contains a Python 3 conversion of the FunSeq2 Perl scripts.

## Files Created

1. **`code/funseq2.py`** - Main script (converted from `funseq2.pl`)
2. **`code/Funseq_SNV.py`** - SNV analysis module (converted from `Funseq_SNV.pm`)
3. **`code/Funseq_Indel.py`** - Indel analysis module (converted from `Funseq_Indel.pm`)
4. **`requirements.txt`** - Python dependencies (minimal - mostly uses standard library)
5. **`CONVERSION_NOTES.md`** - Detailed conversion notes

## Quick Start

### Prerequisites

The Python version still requires the same external command-line tools as the Perl version:

- `bedtools` (for `intersectBed`)
- `tabix`
- `bigWigAverageOverBed` (from UCSC tools)
- `snpMapper` and `indelMapper` (VAT tools)
- `TFMpvalue-sc2pv` (for motif scoring)
- `fastaFromBed` (from BEDTools)
- `R` (for differential expression analysis)

Make sure these are installed and available in your PATH.

### Usage

The Python version uses the same command-line interface as the Perl version:

```bash
python code/funseq2.py <infile> <maf> <genome_mode> <informat> <outformat> \
    <nc_mode> <output_path> <num_per_run> <cancer_type> <score_cut> \
    <weight_mode> <user_anno_dir> <recur_db_use> <sv_length_cut> \
    [gene_list] [--expression EXPR] [--class CLASS] [--exp_format FORMAT] \
    [--config CONFIG_FILE]
```

**Arguments:**
- `infile`: Input variants file (BED or VCF format)
- `maf`: Minor allele frequency threshold (0.0-1.0)
- `genome_mode`: 1 for somatic, 2 for germline/personal genome
- `informat`: Input format (`bed` or `vcf`)
- `outformat`: Output format (`bed` or `vcf`)
- `nc_mode`: 1 to do non-coding only, 0 for both coding and non-coding
- `output_path`: Output directory path
- `num_per_run`: Number of genomes per run (for parallel processing)
- `cancer_type`: Cancer type for recurrence analysis
- `score_cut`: Non-coding candidate score cutoff
- `weight_mode`: 0 for unweighted, 1 for weighted scoring
- `user_anno_dir`: Directory for user-specific annotations
- `recur_db_use`: Recurrence database use flag
- `sv_length_cut`: Structural variant length cutoff (integer or 'inf')
- `gene_list`: (Optional) Gene list file
- `--expression`: (Optional) Expression file
- `--class`: (Optional) Sample class file
- `--exp_format`: (Optional) Expression format (`rpkm` or `raw`)
- `--config`: (Optional) Configuration file path (default: `config.txt`)

### Example

```bash
python code/funseq2.py variants.bed 0.01 1 bed bed 0 output/ 1 all 5 1 \
    user_annotations/ 0 inf --config config.txt
```

## Key Changes from Perl Version

1. **Object-Oriented Structure**: Uses Python classes instead of Perl packages
2. **File Handling**: Uses `with` statements for file operations
3. **System Calls**: Uses `subprocess.run()` instead of backticks
4. **Command-Line Arguments**: Uses `argparse` module
5. **String Operations**: Uses Python string methods and f-strings
6. **Data Structures**: Uses Python dictionaries and lists

## Implementation Status

### Fully Implemented
- Configuration file reading
- Input format validation
- Variant filtering (SNV/Indel separation)
- GERP score calculation
- Conserved region annotation
- HOT region annotation
- Sensitive region annotation
- ENCODE annotation
- User annotation
- Gene linking (promoters, enhancers, introns, UTRs)
- Basic output integration

### Partially Implemented (Placeholders)
- **Motif breaking analysis** (`motif_break()` in Funseq_SNV) - Requires full PFM parsing
- **Motif gain analysis** (`motif_gain()` in both modules) - Requires PWM calculation
- **Recurrence analysis** (`recur()` in Funseq_SNV) - Needs cross-sample analysis
- **Parallel processing** - Currently single-threaded
- **Weighted scoring** - Basic structure in place, needs full implementation

### Notes

- Some complex functions are marked with `pass` and need full implementation
- The motif analysis functions require careful porting of the Perl logic
- Parallel processing support is simplified compared to the Perl version
- All external tool calls are preserved (bedtools, tabix, etc.)

## Testing

Before using in production:

1. Test with small input files
2. Compare outputs with the Perl version
3. Verify all external tools are working
4. Test edge cases (empty files, missing annotations, etc.)

## Differences to Be Aware Of

1. **String Escaping**: Python uses `\t` for tabs (not `\\t` in file writes)
2. **File Paths**: Uses `os.path.join()` and `pathlib.Path` for cross-platform compatibility
3. **Error Handling**: Uses exceptions instead of `die` statements
4. **Regular Expressions**: Uses `re` module instead of Perl regex syntax

## Contributing

If you complete any of the placeholder implementations:

1. Test thoroughly
2. Compare outputs with Perl version
3. Update this README with implementation status
4. Add comments explaining complex logic

## License

Same license as the original FunSeq2 Perl implementation.

