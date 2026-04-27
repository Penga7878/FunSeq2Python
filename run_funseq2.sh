#!/bin/bash
# FunSeq2 Runner Script for VCF files
# This script helps run FunSeq2 analysis on VCF files

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default parameters
INPUT_VCF="sample_input.vcf"
MAF=0.01
GENOME_MODE=1
OUTPUT_PATH="./output"
WEIGHT_MODE=1
SCORE_CUT=0.5
CONFIG_FILE="config.txt"

# Parse command line arguments
if [ $# -eq 0 ]; then
    echo -e "${YELLOW}Usage: $0 [options]${NC}"
    echo ""
    echo "Options:"
    echo "  -i, --input VCF_FILE     Input VCF file (default: sample_input.vcf)"
    echo "  -m, --maf MAF            Minor allele frequency threshold (default: 0.01)"
    echo "  -o, --output PATH       Output directory (default: ./output)"
    echo "  -c, --config FILE       Config file (default: config.txt)"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -i sample_input.vcf -m 0.01 -o ./results"
    exit 0
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_VCF="$2"
            shift 2
            ;;
        -m|--maf)
            MAF="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PATH="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Check if input file exists
if [ ! -f "$INPUT_VCF" ]; then
    echo -e "${RED}Error: Input VCF file '$INPUT_VCF' not found!${NC}"
    exit 1
fi

# Check for required dependencies
echo -e "${YELLOW}Checking dependencies...${NC}"

MISSING_DEPS=()

if ! command -v bedtools &> /dev/null && ! command -v intersectBed &> /dev/null; then
    MISSING_DEPS+=("bedtools")
fi

if ! command -v tabix &> /dev/null; then
    MISSING_DEPS+=("tabix")
fi

if ! command -v bigWigAverageOverBed &> /dev/null; then
    MISSING_DEPS+=("bigWigAverageOverBed (UCSC tools)")
fi

if ! command -v fastaFromBed &> /dev/null; then
    MISSING_DEPS+=("fastaFromBed (bedtools)")
fi

if [ ${#MISSING_DEPS[@]} -gt 0 ]; then
    echo -e "${RED}Missing dependencies:${NC}"
    for dep in "${MISSING_DEPS[@]}"; do
        echo -e "  - $dep"
    done
    echo ""
    echo -e "${YELLOW}Installation instructions:${NC}"
    echo "  macOS: brew install bedtools ucsc-tools htslib"
    echo "  Linux: sudo apt-get install bedtools tabix ucsc-tools"
    echo ""
    echo -e "${YELLOW}Note: Some tools (snpMapper, TFMpvalue-sc2pv) may need to be installed separately.${NC}"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    echo -e "${GREEN}All basic dependencies found!${NC}"
fi

# Check Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: python3 not found!${NC}"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_PATH"

# Run FunSeq2
echo ""
echo -e "${GREEN}Running FunSeq2 analysis...${NC}"
echo "Input: $INPUT_VCF"
echo "Output: $OUTPUT_PATH"
echo ""

python3 code/funseq2.py \
    "$INPUT_VCF" \
    "$MAF" \
    "$GENOME_MODE" \
    vcf \
    vcf \
    0 \
    "$OUTPUT_PATH" \
    1 \
    "general" \
    "$SCORE_CUT" \
    "$WEIGHT_MODE" \
    ./user_annotations \
    "no" \
    inf \
    --config "$CONFIG_FILE"

echo ""
echo -e "${GREEN}Analysis complete!${NC}"
echo "Results are in: $OUTPUT_PATH"
echo "  - Output.vcf: Main results file"
echo "  - Output.indel.vcf: Indel results (if any)"
echo "  - Error.log: Error log file"



