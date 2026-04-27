# Variant Data Structure in FunSeq2

## How Variants Are Stored

**Important:** Variants are NOT stored as objects. Instead, they use a **dictionary-based approach** where:

1. Each variant is identified by a **unique string ID** (variant identifier)
2. Multiple dictionaries store different attributes, all indexed by the same variant ID
3. Think of it like a "distributed object" - the variant's attributes are spread across multiple dictionaries

## Variant ID Format

Each variant has a unique identifier string:

```
"chr1\t100\tG\tT"
```

**Format:** `{chromosome}\t{position}\t{reference_allele}\t{alternate_allele}`

**Example:**
- `"chr1\t100\tA\tT"` = Chromosome 1, position 100, A→T mutation
- `"chr2\t54321\tG\tC"` = Chromosome 2, position 54321, G→C mutation

Note: Sometimes just `chr\tpos` is used (for GERP scores), but full ID includes ref/alt.

## The Dictionary Structure

The `Funseq_SNV` class contains multiple dictionaries, all using variant IDs as keys:

```python
class Funseq_SNV:
    def __init__(self):
        # Each dictionary uses variant_id as key
        self.DES = {}      # variant_id → description string
        self.VCF = {}      # variant_id → VCF line
        self.GERP = {}     # variant_id → GERP score (float)
        self.CONS = {}     # variant_id → 1 (flag: is conserved)
        self.HOT = {}      # variant_id → {hot_region_id: 1, ...}
        self.SEN = {}      # variant_id → 1 (flag: is sensitive)
        self.USEN = {}     # variant_id → 1 (flag: is ultra-sensitive)
        self.ANNO = {}     # variant_id → {annotation_tag: 1, ...}
        self.USER = {}     # variant_id → {user_anno: 1, ...}
        self.GENE = {}     # variant_id → {gene_name: {region_type: 1, ...}, ...}
        self.HUB = {}      # variant_id → {hub_info: 1, ...}
        self.MOTIFBR = {}  # variant_id → motif_breaking_string
        self.MOTIFG = {}   # variant_id → motif_gain_string
        self.VAT = {}      # variant_id → VAT_annotation_string
        self.SELECTION = {} # variant_id → 1 (flag: under selection)
        self.NET_PROB = {} # variant_id → network_probability (float)
        self.NC = {}       # variant_id → 1 (flag: is non-coding)
```

## Visual Representation

Think of it like this - all dictionaries share the same keys (variant IDs):

```
Variant ID: "chr1\t100\tG\tT"
    ↓
    ├─→ self.DES["chr1\t100\tG\tT"] = "chr1\t99\t100\tG\tT"
    ├─→ self.GERP["chr1\t99"] = 2.5
    ├─→ self.CONS["chr1\t100\tG\tT"] = 1  ✓ (in conserved region)
    ├─→ self.ANNO["chr1\t100\tG\tT"] = {"DHS(K562|chr1:99-101)": 1, "TFP(p53|chr1:99-101)": 1}
    ├─→ self.GENE["chr1\t100\tG\tT"] = {"TP53": {"Promoter": 1}, "BRCA1": {"Distal": 1}}
    ├─→ self.MOTIFBR["chr1\t100\tG\tT"] = "p53#chr1:99-101#+"
    ├─→ self.HOT["chr1\t100\tG\tT"] = {"hot_region_1": 1}
    └─→ self.SEN["chr1\t100\tG\tT"] = 1  ✓ (in sensitive region)
```

## Example: Accessing a Variant's Attributes

To get all information about a variant, you'd query multiple dictionaries:

```python
variant_id = "chr1\t100\tG\tT"

# Check if variant exists
if variant_id in self.DES:
    # Get basic description
    description = self.DES[variant_id]  # "chr1\t99\t100\tG\tT"
    
    # Get GERP score (uses shorter ID)
    id_short = "\t".join(variant_id.split("\t")[:2])  # "chr1\t100"
    gerp_score = self.GERP.get(id_short, ".")  # 2.5 or "."
    
    # Check if in conserved region
    is_conserved = variant_id in self.CONS  # True/False
    
    # Get all ENCODE annotations
    annotations = list(self.ANNO.get(variant_id, {}).keys())  # ["DHS(...)", "TFP(...)"]
    
    # Get gene associations
    genes = list(self.GENE.get(variant_id, {}).keys())  # ["TP53", "BRCA1"]
    
    # Get motif breaking info
    motif_break = self.MOTIFBR.get(variant_id, None)  # "p53#..." or None
```

## How Attributes Are Added

Attributes are added incrementally as the pipeline processes variants:

### Step 1: Initial Creation (snv_filter)
```python
variant_id = "chr1\t100\tG\tT"
self.DES[variant_id] = "chr1\t99\t100\tG\tT"  # Basic description
```

### Step 2: GERP Score (gerp_score)
```python
id_short = "chr1\t100"
self.GERP[id_short] = 2.5  # Conservation score
```

### Step 3: Conserved Regions (conserved)
```python
if variant_in_conserved_region:
    self.CONS[variant_id] = 1  # Flag: is conserved
```

### Step 4: ENCODE Annotations (read_encode)
```python
self.ANNO[variant_id] = {}  # Initialize as dict
self.ANNO[variant_id]["DHS(K562|chr1:99-101)"] = 1
self.ANNO[variant_id]["TFP(p53|chr1:99-101)"] = 1
```

### Step 5: Gene Linking (gene_link)
```python
self.GENE[variant_id] = {}  # Initialize as nested dict
self.GENE[variant_id]["TP53"] = {}
self.GENE[variant_id]["TP53"]["Promoter"] = 1
```

## Different Dictionary Value Types

### 1. Simple Flags (Boolean-like)
- `self.CONS[variant_id] = 1`  # Just indicates presence
- `self.SEN[variant_id] = 1`
- `self.USEN[variant_id] = 1`

### 2. Single Values
- `self.GERP[id_short] = 2.5`  # Float score
- `self.NET_PROB[variant_id] = 0.75`  # Probability
- `self.VAT[variant_id] = "nonsynonymous"`  # String annotation

### 3. Dictionaries (Multiple Values)
- `self.ANNO[variant_id] = {"tag1": 1, "tag2": 1}`  # Multiple annotations
- `self.HOT[variant_id] = {"hot1": 1, "hot2": 1}`  # Multiple HOT regions
- `self.HUB[variant_id] = {"hub1": 1, "hub2": 1}`  # Multiple hubs

### 4. Nested Dictionaries
- `self.GENE[variant_id]["TP53"]["Promoter"] = 1`  # Gene → Region type
- `self.USER[variant_id]["CUSTOM(region)"] = 1`  # User annotations

### 5. Strings (Complex Data)
- `self.MOTIFBR[variant_id] = "p53#chr1:99-101#+#1#0.5#0.8"`  # Delimited string
- `self.MOTIFG[variant_id] = "p53#chr1:99-101#+#0.5"`  # Delimited string

## Why This Design?

**Advantages:**
1. **Memory Efficient**: Only stores attributes that exist (sparse storage)
2. **Fast Lookup**: Dictionary access is O(1)
3. **Flexible**: Easy to add new attribute types
4. **Compatible**: Matches the Perl hash-based structure

**Trade-offs:**
- Not object-oriented (would need to create Variant class)
- Attributes scattered across multiple dictionaries
- Need to check multiple dicts to get all attributes

## Conceptual Model

Think of each variant as a "distributed object":

```
Variant "chr1\t100\tG\tT":
┌─────────────────────────────────────┐
│  Basic Info (DES)                   │ → "chr1\t99\t100\tG\tT"
│  Conservation (GERP)                │ → 2.5
│  Location Features:                 │
│    - Conserved (CONS)               │ → ✓
│    - Sensitive (SEN)                │ → ✓
│    - HOT Region (HOT)               │ → {"hot1": 1}
│  Functional Annotations:            │
│    - ENCODE (ANNO)                  │ → {"DHS(...)": 1, "TFP(...)": 1}
│    - Motif Breaking (MOTIFBR)       │ → "p53#..."
│    - Motif Gain (MOTIFG)            │ → "p53#..."
│  Gene Associations:                 │
│    - Genes (GENE)                   │ → {"TP53": {"Promoter": 1}}
│    - Network Hubs (HUB)             │ → {"hub1": 1}
│  Functional Impact:                 │
│    - VAT (for coding)               │ → "nonsynonymous"
│    - Selection (SELECTION)          │ → ✓
└─────────────────────────────────────┘
```

All of this data is stored across multiple dictionaries, but conceptually represents one variant's complete annotation profile.

