# Detailed Explanation of funseq2.py

This document provides a comprehensive explanation of how the FunSeq2 Python script works, breaking down each component and the overall workflow.

## Overview

FunSeq2 is a pipeline for annotating and prioritizing cancer somatic mutations. It takes variant files (BED or VCF format) and enriches them with functional annotations, scores them, and identifies potential driver mutations.

## High-Level Workflow

```
Input Variants → Filter → Annotate → Score → Output
     ↓              ↓         ↓         ↓        ↓
  (BED/VCF)   (MAF filter) (ENCODE,  (Weighted) (BED/VCF
                          GERP, etc.)  scoring   with scores)
```

## Code Structure Breakdown

### 1. Imports and Setup (Lines 1-19)

```python
import os, sys, re, argparse, subprocess, multiprocessing, glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Funseq_SNV import Funseq_SNV
from Funseq_Indel import Funseq_Indel
```

**Purpose:**
- Standard library imports for file operations, regex, command-line parsing, and system calls
- Imports the two main analysis modules: `Funseq_SNV` (for single nucleotide variants) and `Funseq_Indel` (for insertions/deletions)

### 2. Configuration Reading (Lines 22-41)

```python
def read_config(config_file: str) -> Tuple[str, Dict[str, str]]:
```

**What it does:**
- Reads `config.txt` which contains paths to all required data files
- The config file has two types of entries:
  - `file_path=data_context` - Base directory for all data files
  - `tgp_snp=1kg.phase1.snp.bed.gz` - Individual file paths (relative to file_path)

**Example config.txt:**
```
file_path=data_context
tgp_snp=1kg.phase1.snp.bed.gz
gerp_file=All_hg19_RS.bw
encode_annotation=ENCODE.annotation.gz
...
```

**Returns:**
- `file_path`: Base directory path
- `variables`: Dictionary mapping variable names to full file paths
  - Example: `{'tgp_snp': 'data_context/1kg.phase1.snp.bed.gz', ...}`

### 3. Parameter Validation (Lines 44-63)

```python
def validate_parameters(maf, genome_mode, informat, outformat, sv_length_cut, exp_format):
```

**What it checks:**
- **MAF (Minor Allele Frequency)**: Must be between 0.0 and 1.0
- **Genome Mode**: Must be 1 (somatic) or 2 (germline/personal)
- **Input/Output Format**: Must be 'bed' or 'vcf'
- **SV Length Cut**: Must be an integer or 'inf'
- **Expression Format**: If provided, must be 'rpkm' or 'raw'

**Why it matters:** Prevents runtime errors by catching invalid parameters early.

### 4. Main Analysis Pipeline (Lines 72-227)

This is the core function that processes a single input file. Let's break it down step by step:

#### 4.1 Initialization (Lines 80-101)

```python
def main_analysis(infile, output_path, tag, variables, ...):
```

**Setup:**
- Checks input file exists and is not empty
- Cleans up any previous outputs for this sample (`tag`)
- Creates error log file (`{tag}.err`)
- Defines intermediate file paths:
  - `out_snv_filter`: SNVs after MAF filtering
  - `out_nc`: Non-coding variants
  - `out_cds`: Coding variants
  - `out_vat`: VAT (Variant Annotation Tool) output
  - `out_motif`: Motif analysis results
  - `out_indel`: Indels separated from SNVs

#### 4.2 Gene List Filtering (Optional, Lines 103-124)

**If a gene list is provided:**
- Extracts genomic regions (CDS, promoter, intron, UTR, enhancer) for genes in the list
- Filters input variants to only those overlapping these regions
- **Purpose:** Focus analysis on specific genes of interest

**How it works:**
1. Uses `grep` to find gene names in annotation files
2. Extracts chromosome, start, end coordinates
3. Uses `intersectBed` to find variants overlapping these regions

#### 4.3 Format Validation (Lines 130-136)

```python
status = Funseq_SNV.format_check(infile, informat, err_file)
```

**What it checks:**
- **BED format:** Must have at least 5 columns: chr, start, end, ref, alt
- **VCF format:** Must follow VCF 4.0 specification

**Purpose:** Ensures input file is properly formatted before processing.

#### 4.4 Step 1: SNV Filtering (Lines 138-143)

```python
data = Funseq_SNV()
data.snv_filter(infile, informat, variables['tgp_snp'], maf, out_snv_filter, out_indel, sv_length_cut)
```

**What happens:**
1. **Separates SNVs from Indels:**
   - SNVs: Single base changes (A→T, G→C, etc.)
   - Indels: Insertions/deletions (saved to `out_indel`)

2. **Filters against 1000 Genomes Project:**
   - Uses `tabix` to query variants in the 1000 Genomes database
   - Removes variants with MAF ≥ threshold (common variants)
   - **Rationale:** Common variants are less likely to be disease-causing

3. **Output:**
   - `out_snv_filter`: Rare SNVs (potential drivers)
   - `out_indel`: Indels (processed separately later)

**Key concept:** Only rare variants pass this filter. Common variants are filtered out as they're likely benign polymorphisms.

#### 4.5 GERP Score Calculation (Line 152)

```python
data.gerp_score(variables['gerp_file'], out_snv_filter)
```

**What is GERP?**
- GERP++ (Genomic Evolutionary Rate Profiling) measures evolutionary conservation
- Higher scores = more conserved = more likely to be functionally important

**How it works:**
1. Uses `bigWigAverageOverBed` to extract GERP scores for each variant position
2. Stores scores in `data.GERP` dictionary
3. **Purpose:** Variants in highly conserved regions are prioritized

#### 4.6 Step 2: Coding vs Non-Coding Separation (Lines 154-166)

```python
subprocess.run(f"intersectBed -u -a {out_snv_filter} -b {variables['cds']} > {out_cds}")
subprocess.run(f"intersectBed -v -a {out_snv_filter} -b {out_cds} > {out_nc}")
```

**What happens:**
- **Coding variants (`out_cds`):** Overlap with CDS (Coding DNA Sequence) regions
- **Non-coding variants (`out_nc`):** Everything else (promoters, enhancers, introns, intergenic)

**Why separate?**
- Coding variants: Directly affect protein sequence → analyzed with VAT
- Non-coding variants: Affect gene regulation → analyzed with ENCODE annotations

#### 4.7 ENCODE Annotation (Line 164)

```python
data.read_encode(out_nc, variables['encode_annotation'])
```

**What is ENCODE?**
- Encyclopedia of DNA Elements - functional annotations of the genome
- Includes: DHS (DNase Hypersensitive Sites), TFP (Transcription Factor Peaks), Enhancers, etc.

**What it does:**
- Uses `intersectBed` to find which non-coding variants overlap ENCODE annotations
- Stores annotations in `data.ANNO` dictionary
- **Example annotation:** `DHS(K562|chr1:1000-2000)` means variant is in a DNase hypersensitive site in K562 cells

#### 4.8 User Annotations (Line 166)

```python
data.user_annotation(out_snv_filter, user_anno_dir)
```

**Purpose:** Allows users to add custom annotations (e.g., custom regulatory regions, disease-associated regions)

**How it works:**
- Scans `user_anno_dir` for BED files
- Intersects variants with each custom annotation file
- Adds annotations to `data.USER` dictionary

#### 4.9 Step 3: Non-Coding Analysis (Lines 168-180)

This is where non-coding variants get enriched with functional information:

```python
data.conserved(out_snv_filter, variables['conserved'])      # Ultra-conserved regions
data.sensitive(out_nc, variables['sensitive'])              # Sensitive/ultra-sensitive regions
data.hot_region(out_nc, variables['hot_file'])              # Highly occupied target regions
data.motif_break(...)                                        # Transcription factor motif disruption
data.gene_link(...)                                          # Link to nearby genes
data.motif_gain(...)                                         # Novel motif creation
```

**Each annotation type:**

1. **Conserved Regions:**
   - Ultra-conserved elements across species
   - Variants here are highly prioritized

2. **Sensitive Regions:**
   - Genomic regions sensitive to mutations
   - Ultra-sensitive regions get even higher priority

3. **HOT Regions:**
   - Highly Occupied Target regions (many transcription factors bind)
   - Indicates important regulatory regions

4. **Motif Breaking:**
   - Checks if variant disrupts existing transcription factor binding sites
   - Uses Position Frequency Matrices (PFM) to calculate binding affinity
   - **Example:** Variant might break a p53 binding site

5. **Gene Linking:**
   - Links non-coding variants to nearby genes via:
     - Promoters (within 2.5kb upstream)
     - Enhancers (distal regulatory elements)
     - Introns, UTRs
   - Also calculates network centrality (hub genes)

6. **Motif Gain:**
   - Checks if variant creates NEW transcription factor binding sites
   - **Example:** Variant might create a new binding site for a cancer-related TF

#### 4.10 Step 4: Coding Analysis (Lines 182-192)

```python
Funseq_SNV.coding(out_cds, variables['coding_interval'], variables['coding_fasta'], ...)
```

**What is VAT?**
- Variant Annotation Tool (snpMapper/indelMapper)
- Predicts functional impact of coding variants

**What it does:**
1. Converts variants to VCF format
2. Runs `snpMapper` to annotate:
   - Synonymous vs non-synonymous
   - Premature stop codons
   - Frameshifts
   - etc.
3. Also checks:
   - Network hubs (highly connected genes)
   - Negative selection (genes under purifying selection)

**Output:** Variants annotated with functional predictions (e.g., "nonsynonymous", "prematureStop")

#### 4.11 Step 5: Integration and Scoring (Lines 194-196)

```python
data.intergrate(outformat, tag, out_nc, out_vat, ...)
```

**What happens:**
1. **Reads all annotations** collected so far
2. **Calculates scores:**
   - **Coding variants:** Based on functional impact (nonsynonymous +1, premature stop +2, etc.)
   - **Non-coding variants:** Based on weighted sum of features:
     - ENCODE annotation: +weight
     - Sensitive region: +weight
     - Ultra-sensitive: +higher weight
     - Conserved: +weight
     - GERP > 2: +weight
     - HOT region: +weight
     - Motif breaking: +weight (scaled by probability)
     - Motif gain: +weight (scaled by probability)
     - Network hub: +weight (scaled by centrality)
     - Gene association: +weight

3. **Outputs ranked variants:**
   - Coding variants first (sorted by score, highest first)
   - Then non-coding variants (sorted by score, highest first)
   - Format: BED or VCF with all annotations in INFO field

#### 4.12 Indel Analysis (Lines 198-217)

Similar pipeline for indels:
1. GERP scores
2. Annotations (conserved, HOT, sensitive, ENCODE)
3. Gene linking
4. Coding analysis (if applicable)
5. Motif gain analysis
6. Integration and scoring

**Note:** Indels are processed separately because they require different handling (e.g., different motif analysis).

### 5. Main Function (Lines 230-376)

This orchestrates the entire pipeline:

#### 5.1 Argument Parsing (Lines 232-256)

Uses `argparse` to parse command-line arguments. All 14+ required arguments plus optional ones.

#### 5.2 Configuration and Setup (Lines 271-316)

1. Reads config file
2. Auto-detects GENCODE version (finds latest `.promoter.bed` file)
3. Constructs full paths to all annotation files
4. Validates all required files exist

#### 5.3 Output File Setup (Lines 318-328)

Defines output files:
- `Output.{bed|vcf}`: Detailed results
- `Recur.Summary`: Recurrence analysis
- `Candidates.Summary`: Driver candidates
- `Error.log`: Error messages
- `Output.indel.{bed|vcf}`: Indel results

#### 5.4 Differential Expression Analysis (Optional, Lines 330-336)

If expression data is provided:
- Runs R script `differential_gene_expression.r`
- Identifies up/down-regulated genes
- Used later to prioritize variants in differentially expressed genes

#### 5.5 Process Input Files (Lines 338-363)

Currently handles single file case:
1. Extracts sample tag from filename
2. Calls `main_analysis()` for the file
3. Moves/renames output files to final locations

**Note:** Parallel processing for multiple files is simplified (TODO in code).

#### 5.6 Recurrence Analysis (Lines 365-372)

```python
Funseq_SNV.recur(file_detail, file_recur, file_driver, ...)
```

**What is recurrence?**
- Variants that appear in multiple samples are more likely to be drivers
- Can recur at:
  - Same position (same variant in multiple samples)
  - Same gene (different variants in same gene)
  - Same regulatory element (variants in same enhancer/promoter)

**What it does:**
1. Scans all variants across all samples
2. Identifies recurrent elements
3. Adds recurrence information to output
4. Creates `Candidates.Summary` with high-confidence drivers

## Data Flow Diagram

```
Input Variants (BED/VCF)
    ↓
[Format Check]
    ↓
[MAF Filter] → Common variants removed
    ↓
[SNV/Indel Split]
    ↓
┌─────────────────┬─────────────────┐
│   SNVs         │    Indels       │
│                 │                 │
│ [GERP Score]    │ [GERP Score]    │
│ [Coding/NC]     │ [Annotations]   │
│                 │ [Gene Link]     │
│ ┌─────────────┐ │ [Coding]        │
│ │  Coding     │ │ [Motif Gain]   │
│ │  [VAT]      │ │ [Integrate]    │
│ │  [Score]    │ │                 │
│ └─────────────┘ │                 │
│                 │                 │
│ ┌─────────────┐ │                 │
│ │ Non-coding  │ │                 │
│ │ [ENCODE]    │ │                 │
│ │ [Conserved] │ │                 │
│ │ [Sensitive] │ │                 │
│ │ [HOT]       │ │                 │
│ │ [Motif]     │ │                 │
│ │ [Gene Link] │ │                 │
│ │ [Score]     │ │                 │
│ └─────────────┘ │                 │
│                 │                 │
│ [Integrate]     │                 │
└─────────────────┴─────────────────┘
    ↓                    ↓
Output.{format}    Output.indel.{format}
    ↓
[Recurrence Analysis]
    ↓
Candidates.Summary
```

## Key Data Structures

### Funseq_SNV Object Attributes

- `DES`: Variant descriptions (chr, pos, ref, alt)
- `GERP`: GERP++ conservation scores
- `CONS`: Ultra-conserved region flags
- `HOT`: Highly occupied target regions
- `SEN/USEN`: Sensitive/ultra-sensitive regions
- `ANNO`: ENCODE annotations
- `USER`: User-defined annotations
- `GENE`: Gene associations (promoter, enhancer, etc.)
- `HUB`: Network hub information
- `MOTIFBR`: Motif breaking events
- `MOTIFG`: Motif gain events
- `VAT`: Variant annotation tool results
- `SELECTION`: Negative selection flags

## Scoring System

### Coding Variants
- Base score: 0
- +1: Non-synonymous mutation
- +2: Premature stop codon
- +1: Network hub gene
- +1: Under negative selection
- +1: In conserved region
- +1: GERP > 2

### Non-Coding Variants
Uses weighted scoring (if `weight_mode=1`) or unweighted (if `weight_mode=0`):

**Unweighted (simple counting):**
- +1 for each feature present (ENCODE, sensitive, conserved, HOT, motif, etc.)

**Weighted (from weight_file):**
- Each feature has a weight (can be formula like `value * 0.5`)
- More sophisticated scoring based on feature importance

## Output Format

### BED Format
```
chr    start    end    ref    alt    sample    GERP;CDS;VAT;HUB;SELECTION;ENCODE;HOT;MOTIF;SEN;USEN;CONS;GENE;USER;SCORE
```

### VCF Format
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    100  .   A    T    .     .       SAMPLE=sample1;GERP=2.5;CDS=No;NCENC=DHS(...);NCDS=5.2
```

INFO field contains all annotations separated by semicolons.

## Important Concepts

1. **MAF Filtering:** Removes common variants (likely benign)
2. **Conservation:** Variants in conserved regions are prioritized
3. **Regulatory Impact:** Non-coding variants can affect gene expression
4. **Motif Analysis:** Checks both disruption and creation of TF binding sites
5. **Network Analysis:** Prioritizes variants in hub genes (highly connected)
6. **Recurrence:** Variants appearing in multiple samples are more likely drivers
7. **Scoring:** Combines all features into a single score for prioritization

## Limitations and TODOs

1. **Parallel Processing:** Currently simplified (Perl version uses fork manager)
2. **Motif Functions:** `motif_break()` and `motif_gain()` are placeholders
3. **Recurrence Analysis:** `recur()` function is placeholder
4. **Error Handling:** Could be more comprehensive

This pipeline is designed to identify the most functionally significant variants from large cancer sequencing datasets, prioritizing those most likely to be driver mutations.

