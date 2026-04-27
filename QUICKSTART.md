# Quick Start Guide for FunSeq2 with VCF Files

## Setup Complete! ✅

Your data context has been organized in the `data_context/` directory. The following files are set up:

- ✅ 1000 Genomes SNP data
- ✅ GERP scores
- ✅ ENCODE annotations
- ✅ Motif files
- ✅ Gene annotations (GENCODE v19)
- ✅ Regulatory networks
- ✅ Gene lists and cancer recurrence data

## Missing Files

⚠️ **Reference Genome**: The `reference_genome.fa` file is missing. This is a large file (~3GB) needed for motif gain analysis. The analysis will run without it, but motif gain features will be skipped.

To download the reference genome:
```bash
# Option 1: Download from UCSC (hg19/GRCh37)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xzf chromFa.tar.gz
cat chr*.fa > data_context/reference_genome.fa
rm chr*.fa chromFa.tar.gz

# Option 2: Use a pre-indexed version
# Check if you have hg19.fa or similar in your system
```

## Required Dependencies

You need to install the following command-line tools:

### macOS (using Homebrew):
```bash
brew install bedtools ucsc-tools htslib
```

### Linux (Ubuntu/Debian):
```bash
sudo apt-get install bedtools tabix ucsc-tools
```

### Additional Tools (may need separate installation):
- `snpMapper` - VAT tool for coding variant annotation
- `TFMpvalue-sc2pv` - For motif scoring (optional if motif gain is skipped)

## Running the Analysis

### Option 1: Use the Helper Script (Recommended)
```bash
./run_funseq2.sh -i sample_input.vcf -m 0.01 -o ./output
```

### Option 2: Run Directly
```bash
python3 code/funseq2.py \
    sample_input.vcf \
    0.01 \
    1 \
    vcf \
    vcf \
    0 \
    ./output \
    1 \
    "general" \
    0.5 \
    1 \
    ./user_annotations \
    "no" \
    inf \
    --config config.txt
```

## Parameters Explained

- `sample_input.vcf` - Your input VCF file
- `0.01` - MAF threshold (filters common variants)
- `1` - Genome mode (1=somatic, 2=germline)
- `vcf` - Input format
- `vcf` - Output format
- `0` - nc_mode (0=both coding+non-coding, 1=non-coding only)
- `./output` - Output directory
- `0.5` - Score cutoff for non-coding candidates
- `1` - Weight mode (1=weighted, 0=unweighted)

## Output Files

After running, you'll find in the output directory:
- `Output.vcf` - Main results with annotations and scores
- `Output.indel.vcf` - Indel results (if any)
- `Error.log` - Error log
- `Recur.Summary` - Recurrence analysis summary
- `Candidates.Summary` - Candidate driver variants

## Troubleshooting

1. **"intersectBed: command not found"**
   - Install bedtools: `brew install bedtools` (macOS) or `sudo apt-get install bedtools` (Linux)

2. **"tabix: command not found"**
   - Install htslib: `brew install htslib` (macOS) or `sudo apt-get install tabix` (Linux)

3. **"bigWigAverageOverBed: command not found"**
   - Install UCSC tools: `brew install ucsc-tools` (macOS) or download from UCSC

4. **Reference genome warning**
   - This is expected if you haven't downloaded the reference genome yet
   - Analysis will proceed but motif gain analysis will be skipped

## Next Steps

1. Install the required dependencies
2. (Optional) Download the reference genome for full motif gain analysis
3. Run the analysis using the helper script or direct command
4. Check the output files for annotated variants



