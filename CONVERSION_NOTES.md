# FunSeq2 Perl to Python Conversion Notes

This document describes the conversion of FunSeq2 from Perl to Python and important differences to be aware of.

## Overview

The Perl scripts have been converted to Python 3, maintaining the same functionality while using Python idioms and standard library modules.

## File Structure

- `code/funseq2.py` - Main script (converted from `funseq2.pl`)
- `code/Funseq_SNV.py` - SNV analysis module (converted from `Funseq_SNV.pm`)
- `code/Funseq_Indel.py` - Indel analysis module (converted from `Funseq_Indel.pm`)

## Key Conversion Changes

### 1. Data Structures

**Perl:**
```perl
my %hash = ();
$hash{$key} = $value;
```

**Python:**
```python
hash = {}
hash[key] = value
```

### 2. Object-Oriented Structure

Perl uses `bless` for object creation. Python uses classes directly:

**Perl:**
```perl
sub new {
    my $class = shift;
    my $self = {};
    bless($self, $class);
    return $self;
}
```

**Python:**
```python
class Funseq_SNV:
    def __init__(self):
        self.DES = {}
        # ... other attributes
```

### 3. File Operations

**Perl:**
```perl
open(IN, $file) || die;
while(<IN>) {
    # process line
}
close IN;
```

**Python:**
```python
with open(file, 'r') as f:
    for line in f:
        # process line
```

### 4. System Commands

Both Perl and Python use subprocess/system calls extensively. The Python version uses `subprocess.run()`:

**Perl:**
```perl
my $output = `command`;
```

**Python:**
```python
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
output = result.stdout
```

### 5. Regular Expressions

**Perl:**
```perl
if ($string =~ /pattern/) {
    # match
}
```

**Python:**
```python
import re
if re.search(r'pattern', string):
    # match
```

### 6. Command-Line Arguments

**Perl:**
```perl
my $arg = $ARGV[0];
```

**Python:**
```python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('arg')
args = parser.parse_args()
```

### 7. String Operations

**Perl:**
```perl
my $joined = join("\t", @array);
my @split = split(/\t/, $string);
```

**Python:**
```python
joined = "\t".join(array)
split = string.split("\t")
```

## Important Notes

### External Dependencies

The Python version still requires the same external command-line tools:
- `bedtools` (for `intersectBed`)
- `tabix`
- `bigWigAverageOverBed` (from UCSC tools)
- `snpMapper` and `indelMapper` (VAT tools)
- `TFMpvalue-sc2pv` (for motif scoring)
- `fastaFromBed` (from BEDTools)
- `R` (for differential expression analysis)

### Incomplete Implementations

Some complex functions are marked with placeholders (`pass`) and need full implementation:
- `motif_break()` in Funseq_SNV - Requires PFM file parsing and motif scoring
- `motif_gain()` in Funseq_SNV - Requires sequence extraction and PWM calculation
- `motif_gain()` in Funseq_Indel - Similar to SNV version but for indels
- `recur()` in Funseq_SNV - Recurrence analysis across samples

These functions require careful porting of the complex logic from the Perl versions.

### Parallel Processing

The Perl version uses `Parallel::ForkManager` for parallel processing. The Python version currently has a simplified single-threaded implementation. For full parallel processing support, consider using Python's `multiprocessing` module.

### String Escaping

When writing to files, be careful with string escaping. In Python f-strings, use single backslashes for tab (`\t`) and newline (`\n`), not double backslashes.

## Testing Recommendations

1. Test with small input files first
2. Compare outputs between Perl and Python versions
3. Verify all external tools are available and in PATH
4. Test edge cases (empty files, missing annotations, etc.)

## Usage

The Python version maintains the same command-line interface as the Perl version:

```bash
python code/funseq2.py <infile> <maf> <genome_mode> <informat> <outformat> \
    <nc_mode> <output_path> <num_per_run> <cancer_type> <score_cut> \
    <weight_mode> <user_anno_dir> <recur_db_use> <sv_length_cut> \
    [gene_list] [--expression EXPR] [--class CLASS] [--exp_format FORMAT]
```

## Future Improvements

1. Complete implementation of motif analysis functions
2. Add proper parallel processing support
3. Add comprehensive error handling
4. Add unit tests
5. Consider using Python libraries (e.g., `pysam` for BAM/VCF handling, `pandas` for data manipulation)

