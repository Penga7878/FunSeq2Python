"""
Funseq_SNV Python Module
Handles SNV (Single Nucleotide Variant) analysis
"""

import os
import re
import subprocess
import math
from typing import Dict, List, Optional, Tuple


class Funseq_SNV:
    """Class for SNV analysis and annotation."""
    
    def __init__(self):
        """Initialize the object with empty dictionaries for data storage."""
        self.DES = {}  # Variant descriptions
        self.VCF = {}  # VCF format data
        self.GERP = {}  # GERP scores
        self.CONS = {}  # Conserved regions
        self.HOT = {}  # HOT regions
        self.SEN = {}  # Sensitive regions
        self.USEN = {}  # Ultra-sensitive regions
        self.ANNO = {}  # Annotations
        self.USER = {}  # User annotations
        self.GENE = {}  # Gene associations
        self.HUB = {}  # Network hubs
        self.MOTIFBR = {}  # Motif breaking
        self.BR_PROB = {}  # Breaking probability
        self.MOTIFG = {}  # Motif gain
        self.G_PROB = {}  # Gain probability
        self.VAT = {}  # VAT annotations
        self.SELECTION = {}  # Selection data
        self.NET_PROB = {}  # Network probability
        self.NC = {}  # Non-coding variants
    
    @staticmethod
    def format_check(input_file: str, format_type: str, err_file: str) -> int:
        """Check input file format - BED or VCF format."""
        err_bed = """Error: incorrect bed format; please prepare your file as : chromosome\tstart\tend\treference_allele\talternative_allele
* Note: the 6th column is saved for sample names if having multiple samples.

e.g.,
chr1  213941196  213941197	G	T	sample1	  ...
chr1  213942363  213942364	A	C	sample1	  ...
chr1  213943530  213943531	T	A	sample1	  ...
"""
        
        err_vcf = """Error: incorrect vcf format; please prepare your file as :

e.g.,
#CHROM POS    ID        REF  ALT     QUAL FILTER INFO ...
2      4370   rs6057    G    A       29   .      ...   ...
"""
        
        status = 0
        
        if format_type.lower() == 'bed':
            bed_pattern = re.compile(r'\S+\s+\d+\s+\d+\s+\S+\s+\S+')
            with open(input_file, 'r') as f:
                for line in f:
                    if not bed_pattern.match(line.strip()):
                        print(line, end='')
                        status = 1
            if status == 1:
                print(err_bed)
                with open(err_file, 'a') as err:
                    err.write(err_bed)
            return status
        
        if format_type.lower() == 'vcf':
            vcf_pattern = re.compile(r'\S+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+')
            with open(input_file, 'r') as f:
                for line in f:
                    if not line.startswith('#') and not vcf_pattern.match(line.strip()):
                        print(line, end='')
                        status = 1
            if status == 1:
                print(err_vcf)
                with open(err_file, 'a') as err:
                    err.write(err_vcf)
            return status
        
        return status
    
    def snv_filter(self, input_file: str, input_format: str, polymorphism: str,
                   maf_cutoff: float, output: str, out_indel: str, sv_length_cut: str):
        """Filter variants based on minor allele frequency threshold."""
        saw = {}
        tmp_file = f"{output}.tmp"
        
        # Preprocess input file (add/remove 'chr' prefix)
        with open(input_file, 'r') as f_in, open(tmp_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                else:
                    line = re.sub(r'^chr', '', line, flags=re.I)
                    if not line.startswith('#'):
                        line = 'chr' + line
                    f_out.write(line)
        
        # Filter using tabix if MAF cutoff is not 1
        if maf_cutoff != 1:
            cmd = f"tabix {polymorphism} -B {tmp_file} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}} $4 >= {maf_cutoff}'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in result.stdout.strip().split('\n'):
                if line:
                    fields = line.split('\t')
                    if len(fields) >= 2:
                        id_key = f"{fields[0]}\t{fields[1]}"
                        saw[id_key] = 1
        
        with open(output, 'w') as out_f:
            if input_format.lower() == 'vcf':
                with open(tmp_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            fields = line.strip().split('\t')
                            if len(fields) >= 5:
                                ref = fields[3]
                                alt = fields[4]
                                if len(ref) == 1 and len(alt) == 1 and re.match(r'[ATCG]', ref, re.I) and re.match(r'[ATCG]', alt, re.I):
                                    id_key = f"{fields[0]}\t{int(fields[1])-1}"
                                    if id_key not in saw:
                                        full_id = f"{fields[0]}\t{int(fields[1])-1}\t{ref}\t{alt}"
                                        self.VCF[full_id] = '\t'.join(fields[:7])
                                        self.DES[full_id] = f"{fields[0]}\t{int(fields[1])-1}\t{fields[1]}\t{ref}\t{alt}"
                                        out_f.write(f"{self.DES[full_id]}\n")
            
            elif input_format.lower() == 'bed':
                with open(out_indel, 'w') as indel_f:
                    with open(tmp_file, 'r') as f:
                        for line in f:
                            fields = line.strip().split('\t')
                            if len(fields) >= 5:
                                ref = fields[3]
                                alt = fields[4]
                                if len(ref) == 1 and len(alt) == 1 and re.match(r'[ATCG]', ref, re.I) and re.match(r'[ATCG]', alt, re.I):
                                    id_key = f"{fields[0]}\t{fields[1]}"
                                    if id_key not in saw:
                                        full_id = f"{fields[0]}\t{fields[1]}\t{ref}\t{alt}"
                                        self.DES[full_id] = '\t'.join(fields[:5])
                                        out_f.write(f"{self.DES[full_id]}\n")
                                else:
                                    # Indel
                                    if sv_length_cut == 'inf' or max(len(ref), len(alt)) <= int(sv_length_cut):
                                        indel_f.write('\t'.join(fields[:5]) + '\n')
        
        os.remove(tmp_file)
    
    def gerp_score(self, gerp_file: str, infile: str):
        """Obtain GERP++ scores for variants."""
        # Initialize GERP scores
        for variant_id in self.DES.keys():
            id_short = '\t'.join(variant_id.split('\t')[:2])
            self.GERP[id_short] = "."
        
        if os.path.isfile(gerp_file) and os.path.getsize(gerp_file) > 0:
            cmd = f"awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,$2,$2+1,$1\":\"$2}}' {infile} | sort -k 1,1 -k 2,2n | uniq | bigWigAverageOverBed {gerp_file} stdin stdout"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = re.split(r'[:|\s]+', line)
                    if len(parts) >= 5:
                        id_key = f"{parts[0]}\t{parts[1]}"
                        score = float(parts[4])
                        if score != 0:
                            self.GERP[id_key] = score
    
    def conserved(self, input_file: str, conserved_file: str):
        """Identify variants in conserved regions."""
        cmd = f"intersectBed -u -a {input_file} -b {conserved_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                    self.CONS[variant_id] = 1
    
    def hot_region(self, input_file: str, hot_file: str):
        """Identify variants in HOT (highly occupied target) regions."""
        cmd = f"intersectBed -a {input_file} -b {hot_file} -wo"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 9:
                    variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                    if variant_id not in self.HOT:
                        self.HOT[variant_id] = {}
                    self.HOT[variant_id][fields[8]] = 1
    
    def sensitive(self, input_file: str, sensitive_file: str):
        """Identify variants in sensitive/ultra-sensitive regions."""
        # Regular sensitive regions
        cmd = f"intersectBed -u -a {input_file} -b {sensitive_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                    self.SEN[variant_id] = 1
        
        # Ultra-sensitive regions
        cmd = f"grep 'Ultra' {sensitive_file} | intersectBed -u -a {input_file} -b stdin"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                    self.USEN[variant_id] = 1
    
    def read_encode(self, input_file: str, encode_annotation: str):
        """Read ENCODE annotations for non-coding variants."""
        cmd = f"intersectBed -a {input_file} -b {encode_annotation} -wo"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split('\t')
                if len(fields) >= 9:
                    variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                    interval = f"{fields[5]}:{fields[6]}-{fields[7]}"
                    anno_parts = fields[8].split('.')
                    
                    if anno_parts and anno_parts[-1] == 'u':
                        tag = f"{anno_parts[0]}({'.'.join(anno_parts[1:-1])})"
                    else:
                        tag = f"{anno_parts[0]}({'.'.join(anno_parts[1:])}|{interval})"
                    
                    if variant_id not in self.ANNO:
                        self.ANNO[variant_id] = {}
                    self.ANNO[variant_id][tag] = 1
    
    def user_annotation(self, input_file: str, anno_dir: str):
        """Add user-specific annotations."""
        if not os.path.isdir(anno_dir):
            return
        
        files = [f for f in os.listdir(anno_dir) if not f.startswith('.')]
        for filename in files:
            filepath = os.path.join(anno_dir, filename)
            if os.path.isfile(filepath):
                cate = filename.split('.')[0].upper()
                cmd = f"intersectBed -a {input_file} -b {filepath} -wo"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                for line in result.stdout.strip().split('\n'):
                    if line:
                        fields = line.split('\t')
                        if len(fields) >= 9:
                            variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                            interval = f"{fields[5]}:{fields[6]}-{fields[7]}"
                            
                            if variant_id not in self.USER:
                                self.USER[variant_id] = {}
                            
                            if len(fields) > 9:
                                tag = f"{cate}({fields[8]}|{interval})"
                            else:
                                tag = f"{cate}({interval})"
                            self.USER[variant_id][tag] = 1
    
    def motif_break(self, input_file: str, ancestral_file: str, bound_motif: str,
                    mode: int, pfm_file: str):
        """Check whether variants break bound motifs."""
        import sys
        import tempfile
        
        # Read in PFM file
        motif = {}
        prev_name = None
        A, C, G, T = 1, 2, 3, 4
        
        try:
            with open(pfm_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        prev_name = line.split('>')[1].split()[0] if '>' in line else line.split()[0]
                    else:
                        info = line.split()
                        if len(info) >= 5:
                            if prev_name not in motif:
                                motif[prev_name] = [{'A': float(info[A]), 'T': float(info[T]), 
                                                      'C': float(info[C]), 'G': float(info[G])}]
                            else:
                                motif[prev_name].append({'A': float(info[A]), 'T': float(info[T]), 
                                                          'C': float(info[C]), 'G': float(info[G])})
        except Exception as e:
            print(f"Error reading PFM file {pfm_file}: {e}", file=sys.stderr)
            return
        
        # Intersect variants with bound motifs
        cmd = (f"intersectBed -a {bound_motif} -b {input_file} -wo | "
               f"sort -k 1,1 -k 2,2n | uniq | "
               f"awk 'BEGIN{{OFS=\"\\t\"}}{{print $8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6,$7}}'")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            input_file_new = result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error running intersectBed: {e}", file=sys.stderr)
            return
        
        # For mode 2 (germline), get ancestral sequences
        ancestral = {}
        if mode == 2:
            # Create temp file for awk processing
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as tmp_bed:
                tmp_bed.write(input_file_new)
                tmp_bed_path = tmp_bed.name
            
            try:
                awk_cmd = ("awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{gsub(\"chr\",\"\",$1); print $1,$2,$3}'")
                fasta_cmd = f"cat {tmp_bed_path} | {awk_cmd} | fastaFromBed -fi {ancestral_file} -bed stdin -fo stdout"
                
                result = subprocess.run(fasta_cmd, shell=True, capture_output=True, text=True, check=True)
                lines = result.stdout.split('\n')
                
                chr_val = None
                pos_val = None
                for line in lines:
                    if line.startswith('>'):
                        parts = re.split(r'[>:\-]', line)
                        if len(parts) >= 3:
                            chr_val = parts[1]
                            pos_val = parts[2]
                            id_short = f"chr{chr_val}\t{pos_val}"
                    elif chr_val and pos_val and line.strip():
                        ancestral[id_short] = line.strip().upper()
            except subprocess.CalledProcessError as e:
                print(f"Warning: Could not get ancestral sequences: {e}", file=sys.stderr)
            finally:
                try:
                    os.unlink(tmp_bed_path)
                except:
                    pass
        
        # Process each intersecting variant
        for input_line in input_file_new.strip().split('\n'):
            if not input_line.strip():
                continue
            
            info = input_line.split('\t')
            if len(info) < 12:
                continue
            
            ref = info[3].upper()
            alt = info[4].upper()
            id_short = f"{info[0]}\t{info[1]}"
            variant_id = f"{info[0]}\t{info[1]}\t{info[3]}\t{info[4]}"
            
            # Determine ancestral and derived alleles
            if mode == 1:  # somatic mode
                AA = ref
                der_al = alt
            else:  # germline mode
                id_key = id_short.replace('chr', '')
                if id_key in ancestral:
                    AA = ancestral[id_key]
                    if AA == ref:
                        der_al = alt
                    elif AA == alt:
                        der_al = ref
                    else:
                        continue
                else:
                    continue
            
            motif_len = int(info[7]) - int(info[6])
            factor = info[11]
            strand = info[10]
            
            # Calculate position in motif
            if strand == "+":
                pos_in_motif = int(info[1]) - int(info[6]) + 1
            elif strand == "-":
                pos_in_motif = int(info[7]) - int(info[1])
                # Reverse complement
                trans_table = str.maketrans('ACGTacgt', 'TGCAtgca')
                AA = AA.translate(trans_table)
                der_al = der_al.translate(trans_table)
            else:
                continue
            
            # Extract motif name (remove _Nmer suffix)
            motif_name = re.split(r'_\d+mer', info[8])[0]
            
            if motif_name not in motif:
                continue
            
            if motif_len != len(motif[motif_name]):
                continue
            
            # Get frequencies
            if pos_in_motif < 1 or pos_in_motif > len(motif[motif_name]):
                continue
            
            derived_allele_freq = motif[motif_name][pos_in_motif - 1].get(der_al, 0.0)
            AA_freq = motif[motif_name][pos_in_motif - 1].get(AA, 0.0)
            
            interval = f"{info[5]}:{info[6]}-{info[7]}"
            
            # Add annotation
            if variant_id not in self.ANNO:
                self.ANNO[variant_id] = {}
            self.ANNO[variant_id][f"TFM({factor}|{info[8]}|{interval})"] = 1
            
            # Check if motif is broken (derived allele frequency < ancestral frequency)
            if derived_allele_freq < AA_freq:
                diff_freq = AA_freq - derived_allele_freq
                
                motif_br_entry = (f"{factor}#{info[8]}#{info[6]}#{info[7]}#{strand}#"
                                 f"{pos_in_motif}#{derived_allele_freq}#{AA_freq}")
                
                if variant_id in self.MOTIFBR:
                    self.MOTIFBR[variant_id] = f"{self.MOTIFBR[variant_id]},{motif_br_entry}"
                    self.BR_PROB[variant_id] = max(diff_freq, self.BR_PROB.get(variant_id, 0))
                else:
                    self.MOTIFBR[variant_id] = motif_br_entry
                    self.BR_PROB[variant_id] = diff_freq
    
    def motif_gain(self, pfm_file: str, reference_file: str, score_file: str,
                   p_cut: float, out_tmp: str):
        """Check whether variants create novel motifs."""
        import sys
        import os
        
        # Read in PFM and transform it to PWM using 20 sequences
        motif = {}
        prev_name = None
        A, C, G, T = 1, 2, 3, 4
        
        try:
            with open(pfm_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        prev_name = line.split('>')[1].split()[0] if '>' in line else line.split()[0]
                    else:
                        info = line.split()
                        if len(info) >= 5:
                            # Convert PFM to PWM: log((count*20+0.25)/21) - log(0.25)
                            pwm_a = math.log((float(info[A]) * 20 + 0.25) / 21) - math.log(0.25)
                            pwm_t = math.log((float(info[T]) * 20 + 0.25) / 21) - math.log(0.25)
                            pwm_c = math.log((float(info[C]) * 20 + 0.25) / 21) - math.log(0.25)
                            pwm_g = math.log((float(info[G]) * 20 + 0.25) / 21) - math.log(0.25)
                            
                            if prev_name not in motif:
                                motif[prev_name] = [{'A': pwm_a, 'T': pwm_t, 'C': pwm_c, 'G': pwm_g}]
                            else:
                                motif[prev_name].append({'A': pwm_a, 'T': pwm_t, 'C': pwm_c, 'G': pwm_g})
        except Exception as e:
            print(f"Error reading PFM file {pfm_file}: {e}", file=sys.stderr)
            return
        
        # Read in score file
        score_lower = {}
        score = {}
        score_upper = {}
        
        try:
            with open(score_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        name = parts[0]
                        cut_off = float(parts[1])
                        p = float(parts[2])
                        
                        if p < p_cut:
                            score_upper[name] = cut_off
                        elif p == p_cut:
                            score[name] = cut_off
                        elif p > p_cut:
                            score_lower[name] = cut_off
        except Exception as e:
            print(f"Error reading score file {score_file}: {e}", file=sys.stderr)
            return
        
        # Retrieve variants in promoter/distal regions and prepare sequences
        variant_ids = []
        ref_alleles = []
        alt_alleles = []
        start_positions = []
        snp_file_content = ""
        
        for variant_id in sorted(self.GENE.keys()):
            switch = False
            for gene in self.GENE[variant_id].keys():
                gene_info = self.GENE[variant_id][gene]
                if "Promoter" in gene_info or "Distal" in gene_info:
                    switch = True
                    break
            
            if switch:
                des_parts = self.DES[variant_id].split()
                if len(des_parts) >= 5:
                    variant_ids.append(variant_id)
                    chr_val = des_parts[0].replace('chr', '')
                    start = int(des_parts[1])
                    ref_alleles.append(des_parts[3].upper())
                    alt_alleles.append(des_parts[4].upper())
                    start_positions.append(start)
                    snp_file_content += f"{chr_val}\t{start-29}\t{start+30}\n"
        
        # Write temp file and extract sequences
        try:
            with open(out_tmp, 'w') as f:
                f.write(snp_file_content)
        except Exception as e:
            print(f"Error writing temp file {out_tmp}: {e}", file=sys.stderr)
            return
        
        if os.path.getsize(out_tmp) == 0:
            return
        
        # Extract sequences using fastaFromBed
        try:
            cmd = f"fastaFromBed -fi {reference_file} -bed {out_tmp} -fo stdout"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            fasta_lines = result.stdout.split('\n')
        except subprocess.CalledProcessError as e:
            print(f"Error extracting sequences: {e}", file=sys.stderr)
            return
        
        if len(fasta_lines) == 0:
            return
        
        # Process each variant
        for idx in range(len(variant_ids)):
            variant_id = variant_ids[idx]
            ref = ref_alleles[idx]
            alt = alt_alleles[idx]
            start = start_positions[idx]
            
            if idx * 2 + 1 >= len(fasta_lines):
                continue
            
            # Helper function for output (defined inside loop to access variant_id, start)
            def out_print(name, i, motif_length, alt_pwm_score, ref_pwm_score, strand):
                """Print motif gain result."""
                diff_score = alt_pwm_score - ref_pwm_score
                
                if strand == "+":
                    tmp = (f"{name}#{start-29+i}#{start+i+motif_length-29}#{strand}#"
                          f"{30-i}#{alt_pwm_score:.3f}#{ref_pwm_score:.3f}")
                else:
                    tmp = (f"{name}#{start-motif_length-i+30}#{start-i+30}#{strand}#"
                          f"{30-i}#{alt_pwm_score:.3f}#{ref_pwm_score:.3f}")
                
                if variant_id in self.MOTIFG:
                    self.MOTIFG[variant_id] = f"{self.MOTIFG[variant_id]},{tmp}"
                    self.G_PROB[variant_id] = max(diff_score, self.G_PROB.get(variant_id, 0))
                else:
                    self.MOTIFG[variant_id] = tmp
                    self.G_PROB[variant_id] = diff_score
            
            # Helper function for sequence scanning (defined inside loop to access ref, alt, variant_id, start)
            def seq_scan(seq_array, strand, ref_allele, alt_allele):
                """Scan sequence for motif binding sites."""
                seq = seq_array
                
                for motif_name in sorted(motif.keys()):
                    motif_length = len(motif[motif_name])
                    
                    # Create PWM directory and files if needed
                    pwm_dir = "pwm"
                    if not os.path.exists(pwm_dir):
                        os.makedirs(pwm_dir)
                    
                    pwm_file_path = os.path.join(pwm_dir, motif_name)
                    
                    # Write PWM file if it doesn't exist
                    if not os.path.exists(pwm_file_path):
                        pwm_lines = []
                        for base in ['A', 'C', 'G', 'T']:
                            pwm_line = []
                            for i in range(motif_length):
                                pwm_line.append(str(motif[motif_name][i][base]))
                            pwm_lines.append('\t'.join(pwm_line))
                        
                        with open(pwm_file_path, 'w') as f:
                            f.write('\n'.join(pwm_lines) + '\n')
                    
                    # Sequence scanning
                    for i in range(30 - motif_length, 30):
                        alt_pwm_score = 0.0
                        ref_pwm_score = 0.0
                        
                        # Calculate alt PWM score
                        for j in range(motif_length):
                            if i + j < len(seq):
                                alt_pwm_score += motif[motif_name][j].get(seq[i + j], 0.0)
                        
                        # Calculate ref PWM score (alt - old base + ref base)
                        ref_pwm_score = (alt_pwm_score - 
                                       motif[motif_name][29 - i].get(alt_allele, 0.0) + 
                                       motif[motif_name][29 - i].get(ref_allele, 0.0))
                        
                        # Check if alt score > ref score (motif gain)
                        if alt_pwm_score > ref_pwm_score:
                            should_print = False
                            
                            if motif_name in score:
                                if alt_pwm_score > score[motif_name] and ref_pwm_score < score[motif_name]:
                                    should_print = True
                            elif ((motif_name in score_upper and ref_pwm_score >= score_upper[motif_name]) or
                                  (motif_name in score_lower and alt_pwm_score <= score_lower[motif_name])):
                                should_print = False
                            elif (motif_name in score_upper and motif_name in score_lower and
                                  alt_pwm_score >= score_upper[motif_name] and 
                                  ref_pwm_score <= score_lower[motif_name]):
                                should_print = True
                            else:
                                # Need to calculate p-values
                                if motif_name in score_upper and alt_pwm_score > score_upper[motif_name]:
                                    # Calculate ref p-value
                                    cmd = (f"TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 "
                                          f"-s {ref_pwm_score} -m {pwm_file_path} -w")
                                    try:
                                        result = subprocess.run(cmd, shell=True, capture_output=True, 
                                                              text=True, check=True)
                                        ref_p = float(result.stdout.split()[2])
                                        if ref_p > p_cut:
                                            should_print = True
                                    except:
                                        pass
                                else:
                                    # Calculate both p-values
                                    cmd_alt = (f"TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 "
                                              f"-s {alt_pwm_score} -m {pwm_file_path} -w")
                                    cmd_ref = (f"TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 "
                                              f"-s {ref_pwm_score} -m {pwm_file_path} -w")
                                    try:
                                        result_alt = subprocess.run(cmd_alt, shell=True, capture_output=True,
                                                                   text=True, check=True)
                                        result_ref = subprocess.run(cmd_ref, shell=True, capture_output=True,
                                                                   text=True, check=True)
                                        alt_p = float(result_alt.stdout.split()[2])
                                        ref_p = float(result_ref.stdout.split()[2])
                                        if alt_p < p_cut and ref_p > p_cut:
                                            should_print = True
                                    except:
                                        pass
                            
                            if should_print:
                                out_print(motif_name, i, motif_length, alt_pwm_score, ref_pwm_score, strand)
            
            extract_seq = fasta_lines[idx * 2 + 1].upper()
            extract_seq_list = list(extract_seq)
            extract_seq_list[29] = alt
            
            # Positive strand
            seq_scan(extract_seq_list, "+", ref, alt)
            
            # Negative strand (reverse complement)
            trans_table = str.maketrans('ATCGatcg', 'TAGCTAGC')
            extract_seq_rc = extract_seq.translate(trans_table)
            extract_seq_list_rc = list(reversed(list(extract_seq_rc)))
            ref_rc = ref.translate(trans_table)
            alt_rc = alt.translate(trans_table)
            extract_seq_list_rc[29] = alt_rc
            
            # Negative strand
            seq_scan(extract_seq_list_rc, "-", ref_rc, alt_rc)
    
    def gene_link(self, input_file: str, promoter: str, distal: str, intron: str,
                  utr: str, network_dir: str):
        """Link non-coding mutations with genes (promoters & regulatory elements)."""
        # Intersect with different gene regions
        cmd_pro = f"intersectBed -a {input_file} -b {promoter} -wo | sort -k 1,1 -k 2,2n | uniq"
        cmd_dis = f"intersectBed -a {input_file} -b {distal} -wo | sort -k 1,1 -k 2,2n | uniq"
        cmd_intron = f"intersectBed -a {input_file} -b {intron} -wo | sort -k 1,1 -k 2,2n | uniq"
        cmd_utr = f"intersectBed -a {input_file} -b {utr} -wo | sort -k 1,1 -k 2,2n | uniq"
        
        # Process each region type
        for cmd, tag in [(cmd_pro, "Promoter"), (cmd_dis, "Distal"), 
                         (cmd_intron, "Intron"), (cmd_utr, "UTR")]:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in result.stdout.strip().split('\n'):
                if line:
                    fields = line.split()
                    if len(fields) >= 9:
                        variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                        gene = fields[8]
                        
                        if variant_id not in self.GENE:
                            self.GENE[variant_id] = {}
                        if gene not in self.GENE[variant_id]:
                            self.GENE[variant_id][gene] = {}
                        self.GENE[variant_id][gene][tag] = 1
        
        # Network analysis (simplified - full implementation would load network files)
        # This would calculate network hub probabilities
    
    @staticmethod
    def coding(snp_input: str, file_interval: str, file_fasta: str, network: str,
               file_selection: str, out: str, out_path: str, nc_mode: int, cds: str) -> int:
        """Coding pipe with VAT (variant annotation tool) analysis."""
        if nc_mode == 0:
            # Convert to VCF format and run snpMapper
            cmd = f"awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print \"##fileformat=VCFv4.0\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\"}}{{print $1,$3,\".\",$4,$5,\".\",\"PASS\",\".\"}}' {snp_input} | snpMapper {file_interval} {file_fasta} | grep '^[^#]'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            vat_out = result.stdout.strip().split('\n')
            
            if not vat_out or not vat_out[0]:
                # Fallback to intersectBed
                cmd = f"intersectBed -a {snp_input} -b {cds} -wo | cut -f 1,3,4,5,9 | sort | uniq"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                with open(out, 'w') as f:
                    for line in result.stdout.strip().split('\n'):
                        if line:
                            fields = line.split()
                            if len(fields) >= 5:
                                joined = '\t'.join(fields[:4])
                                f.write(f"{joined}\t.\t{fields[4]}\t.\n")
                return 1
            else:
                # Process VAT output
                with open(out, 'w') as f:
                    for line in vat_out:
                        if line and not line.startswith('#'):
                            fields = line.split()
                            if len(fields) >= 8:
                                info = fields[7].split(',')[0] if fields[7] else ""
                                if 'VA=' in info:
                                    match = re.search(r'VA=\\d+:(.*?):', info)
                                    if match:
                                        gene = match.group(1)
                                        f.write(f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}\t{fields[7]}\t{gene}\t.\n")
                return 0
        return 0
    
    def intergrate(self, output_format: str, sample: str, nc_snp: str, coding_info: str,
                   gene_dir: str, reg_net: str, out: str, de_data: str, weight_mode: int,
                   weight_file: str):
        """Integrate coding & non-coding outputs with full scoring."""
        nc_score = {}
        cds_score = {}
        gene_info = {}
        
        # Read gene information
        if os.path.isdir(gene_dir):
            files = [f for f in os.listdir(gene_dir) if not f.startswith('.')]
            for filename in files:
                filepath = os.path.join(gene_dir, filename)
                if os.path.isfile(filepath):
                    gene_cate = f"[{filename.split('.')[0]}]"
                    with open(filepath, 'r') as f:
                        for line in f:
                            gene = line.split()[0]
                            if gene not in gene_info:
                                gene_info[gene] = {}
                            gene_info[gene][gene_cate] = 1
        
        # Read differential expression data
        if os.path.isfile(de_data) and os.path.getsize(de_data) > 0:
            with open(de_data, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        gene_list_str = fields[0]
                        cls = fields[1]
                        for gene in gene_list_str.split('|'):
                            if 'benign' in cls.lower():
                                gene_info.setdefault(gene, {})["[down_regulated]"] = 1
                            else:
                                gene_info.setdefault(gene, {})["[up_regulated]"] = 1
        
        # Process regulatory network for cancer gene associations
        cancer_reg = {}
        if os.path.isfile(reg_net):
            with open(reg_net, 'r') as f:
                for line in f:
                    fields = line.split()
                    if len(fields) >= 2:
                        tf = fields[0]
                        gene = fields[1]
                        if gene in gene_info and "[cancer]" in gene_info[gene]:
                            if tf not in cancer_reg:
                                cancer_reg[tf] = {}
                            cancer_reg[tf][gene] = 1
        
        # Add TF regulating cancer gene info
        for tf in cancer_reg:
            info = f"[TF_regulating_known_cancer_gene:{','.join(sorted(cancer_reg[tf].keys()))}]"
            gene_info.setdefault(tf, {})[info] = 1
        
        # Weight function (defined first so it can be used in scoring)
        def weight(file):
            weights = {}
            if os.path.isfile(file):
                with open(file, 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        if len(fields) >= 2:
                            cate = fields[0]
                            func = fields[1]
                            weights[cate] = func
            return weights
        
        # Read non-coding variants
        def read_nc(file):
            if os.path.isfile(file) and os.path.getsize(file) > 0:
                with open(file, 'r') as f:
                    for line in f:
                        fields = line.strip().split()
                        if len(fields) >= 5:
                            variant_id = f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                            self.NC[variant_id] = 1
        
        # Read coding variants and calculate scores
        def read_cds(file):
            if os.path.isfile(file) and os.path.getsize(file) > 0:
                with open(file, 'r') as f:
                    for line in f:
                        line = line.rstrip('\r\n')  # Handle different line endings
                        fields = line.split('\t')
                        if len(fields) >= 4:
                            variant_id = f"{fields[0]}\t{int(fields[1])-1}\t{fields[2]}\t{fields[3]}"
                            vat = fields[4].split(';')[-1] if len(fields) > 4 else "."
                            score = 0
                            self.VAT[variant_id] = vat
                            if len(fields) > 5:
                                self.GENE[variant_id] = fields[5]
                            
                            if 'nonsynonymous' in vat or 'prematureStop' in vat:
                                score += 1
                                if 'prematureStop' in vat:
                                    score += 1
                                if len(fields) > 6 and fields[6] != ".":
                                    self.HUB[variant_id] = fields[6]
                                    score += 1
                                if 'NegativeSelection' in line:
                                    self.SELECTION[variant_id] = 1
                                    score += 1
                                if variant_id in self.CONS:
                                    score += 1
                                
                                id_short = '\t'.join(variant_id.split('\t')[:2])
                                if id_short in self.GERP and self.GERP[id_short] != ".":
                                    try:
                                        gerp_val = float(self.GERP[id_short])
                                        if gerp_val > 2:
                                            score += 1
                                    except (ValueError, TypeError):
                                        pass
                            
                            if score not in cds_score:
                                cds_score[score] = {}
                            cds_score[score][variant_id] = 1
        
        # Weighted scoring scheme
        def nc_score_scheme():
            weights = weight(weight_file)
            for variant_id in self.NC.keys():
                score = 0.0
                id_short = '\t'.join(variant_id.split('\t')[:2])
                
                # ENCODE annotation (only if not SEN, MOTIFBR, or HOT)
                if (variant_id in self.ANNO and 
                    variant_id not in self.SEN and
                    variant_id not in self.MOTIFBR and
                    variant_id not in self.HOT):
                    score = float(weights.get('ANNO', 0))
                
                # Ultra-sensitive or sensitive
                if variant_id in self.USEN:
                    score = float(weights.get('USEN', 0))
                elif variant_id in self.SEN:
                    score = float(weights.get('SEN', 0))
                
                # Motif breaking
                if variant_id in self.MOTIFBR:
                    tmp_score = weights.get('MOTIFBR', '0')
                    if variant_id in self.BR_PROB:
                        tmp_score = tmp_score.replace('value', str(self.BR_PROB[variant_id]))
                    try:
                        score += eval(tmp_score)
                    except:
                        pass
                
                # HOT region
                if variant_id in self.HOT:
                    score += float(weights.get('HOT', 0))
                
                # Conserved or GERP
                if variant_id in self.CONS:
                    score += float(weights.get('CONS', 0))
                elif id_short in self.GERP and self.GERP[id_short] != ".":
                    tmp_score = weights.get('GERP', '0')
                    try:
                        gerp_val = float(self.GERP[id_short])
                        tmp_score = tmp_score.replace('value', str(gerp_val))
                        score += eval(tmp_score)
                    except:
                        pass
                
                # Gene association (only if not HUB or MOTIFG)
                if (variant_id in self.GENE and
                    variant_id not in self.HUB and
                    variant_id not in self.MOTIFG):
                    score += float(weights.get('GENE', 0))
                
                # Network hub
                if variant_id in self.HUB:
                    tmp_score = weights.get('HUB', '0')
                    if variant_id in self.NET_PROB:
                        tmp_score = tmp_score.replace('value', str(self.NET_PROB[variant_id]))
                    try:
                        score += eval(tmp_score)
                    except:
                        pass
                
                # Motif gain
                if variant_id in self.MOTIFG:
                    if variant_id in self.G_PROB:
                        g_prob_val = self.G_PROB[variant_id]
                        import math
                        threshold = math.log(20.25) - math.log(0.25)
                        if g_prob_val < threshold:
                            tmp_score = weights.get('MOTIFG', '0')
                            tmp_score = tmp_score.replace('value', str(g_prob_val))
                            try:
                                score += eval(tmp_score)
                            except:
                                pass
                        else:
                            score += 0.99629183  # Hard-coded value from Perl
                
                if score not in nc_score:
                    nc_score[score] = {}
                nc_score[score][variant_id] = 1
        
        # Unweighted scoring scheme
        def nc_score_scheme_uw():
            for variant_id in self.NC.keys():
                score = 0
                id_short = '\t'.join(variant_id.split('\t')[:2])
                
                if variant_id in self.ANNO:
                    score += 1
                if variant_id in self.SEN:
                    score += 1
                if variant_id in self.USEN:
                    score += 1
                if variant_id in self.CONS:
                    score += 1
                if id_short in self.GERP and self.GERP[id_short] != ".":
                    try:
                        if float(self.GERP[id_short]) > 2:
                            score += 1
                    except:
                        pass
                if variant_id in self.HOT:
                    score += 1
                if variant_id in self.HUB:
                    if variant_id in self.NET_PROB:
                        try:
                            if float(self.NET_PROB[variant_id]) > 0.75:
                                score += 1
                        except:
                            pass
                if variant_id in self.MOTIFBR:
                    score += 1
                if variant_id in self.MOTIFG:
                    score += 1
                if variant_id in self.GENE:
                    score += 1
                
                if score not in nc_score:
                    nc_score[score] = {}
                nc_score[score][variant_id] = 1
        
        # Read data
        read_nc(nc_snp)
        read_cds(coding_info)
        
        # Calculate scores
        if weight_mode == 1:
            nc_score_scheme()
        else:
            nc_score_scheme_uw()
        
        # Output
        def print_bed_format():
            # Coding variants
            for score in sorted(cds_score.keys(), reverse=True):
                for variant_id in sorted(cds_score[score].keys()):
                    if variant_id in self.DES:
                        f_out.write(f"{self.DES[variant_id]}\t{sample}\t")
                        id_short = '\t'.join(variant_id.split('\t')[:2])
                        f_out.write(f"{self.GERP.get(id_short, '.')};")
                        f_out.write("Yes;")
                        f_out.write(f"{self.VAT.get(variant_id, '.')};")
                        
                        if variant_id in self.HUB:
                            hub_val = self.HUB[variant_id]
                            if isinstance(hub_val, dict):
                                f_out.write(f"{','.join(sorted(hub_val.keys()))};")
                            else:
                                f_out.write(f"{hub_val};")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.SELECTION:
                            f_out.write("Yes;")
                        else:
                            f_out.write(".;")
                        
                        f_out.write(".;" * 5)  # 5 empty fields for non-coding features
                        
                        if variant_id in self.CONS:
                            f_out.write("Yes;")
                        else:
                            f_out.write(".;")
                        
                        gene = self.GENE.get(variant_id, ".")
                        if gene in gene_info:
                            f_out.write(f"{gene}{''.join(sorted(gene_info[gene].keys()))};")
                        else:
                            f_out.write(f"{gene};")
                        
                        if hasattr(self, 'USER') and variant_id in self.USER:
                            f_out.write(f"{','.join(sorted(self.USER[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        f_out.write(f"{score};.\n")
            
            # Non-coding variants
            for score in sorted(nc_score.keys(), reverse=True):
                for variant_id in sorted(nc_score[score].keys()):
                    if variant_id in self.DES:
                        f_out.write(f"{self.DES[variant_id]}\t{sample}\t")
                        id_short = '\t'.join(variant_id.split('\t')[:2])
                        f_out.write(f"{self.GERP.get(id_short, '.')};")
                        f_out.write("No;.;")
                        
                        if variant_id in self.HUB:
                            f_out.write(f"{','.join(sorted(self.HUB[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        f_out.write(".;")  # Negative selection
                        
                        if variant_id in self.ANNO:
                            f_out.write(f"{','.join(sorted(self.ANNO[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.HOT:
                            f_out.write(f"{','.join(sorted(self.HOT[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.MOTIFBR or variant_id in self.MOTIFG:
                            if variant_id in self.MOTIFBR:
                                f_out.write(f"MOTIFBR={self.MOTIFBR[variant_id]}")
                                if variant_id in self.MOTIFG:
                                    f_out.write(f",MOTIFG={self.MOTIFG[variant_id]};")
                                else:
                                    f_out.write(";")
                            else:
                                f_out.write(f"MOTIFG={self.MOTIFG[variant_id]};")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.SEN:
                            f_out.write("Yes;")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.USEN:
                            f_out.write("Yes;")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.CONS:
                            f_out.write("Yes;")
                        else:
                            f_out.write(".;")
                        
                        if variant_id in self.GENE:
                            gene_info_str = ""
                            for gene in sorted(self.GENE[variant_id].keys()):
                                info = f"{gene}({'&'.join(sorted(self.GENE[variant_id][gene].keys()))})"
                                if gene in gene_info:
                                    info += ''.join(sorted(gene_info[gene].keys()))
                                gene_info_str += f"{info},"
                            gene_info_str = gene_info_str.rstrip(',')
                            f_out.write(f"{gene_info_str};")
                        else:
                            f_out.write(".;")
                        
                        if hasattr(self, 'USER') and variant_id in self.USER:
                            f_out.write(f"{','.join(sorted(self.USER[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        f_out.write(f".;{score}\n")
        
        def print_vcf_format():
            # Coding variants
            for score in sorted(cds_score.keys(), reverse=True):
                for variant_id in sorted(cds_score[score].keys()):
                    if variant_id in self.VCF:
                        f_out.write(f"{self.VCF[variant_id]}\t")
                    elif variant_id in self.DES:
                        fields = self.DES[variant_id].split('\t')
                        if len(fields) >= 5:
                            f_out.write(f"{fields[0]}\t{fields[2]}\t.\t{fields[3]}\t{fields[4]}\t.\t.\t")
                    else:
                        continue
                    
                    f_out.write(f"SAMPLE={sample};")
                    id_short = '\t'.join(variant_id.split('\t')[:2])
                    f_out.write(f"GERP={self.GERP.get(id_short, '.')};")
                    f_out.write("CDS=Yes;")
                    f_out.write(f"{self.VAT.get(variant_id, '.')};")
                    
                    if variant_id in self.HUB:
                        hub_val = self.HUB[variant_id]
                        if isinstance(hub_val, dict):
                            f_out.write(f"HUB={','.join(sorted(hub_val.keys()))};")
                        else:
                            f_out.write(f"HUB={hub_val};")
                    if variant_id in self.SELECTION:
                        f_out.write("GNEG=Yes;")
                    if variant_id in self.CONS:
                        f_out.write("UCONS=Yes;")
                    
                    gene = self.GENE.get(variant_id, ".")
                    f_out.write(f"GENE={gene};")
                    
                    if gene in gene_info:
                        f_out.write(f"CANG={''.join(sorted(gene_info[gene].keys()))};")
                    
                    if hasattr(self, 'USER') and variant_id in self.USER:
                        f_out.write(f"USER_ANNO={','.join(sorted(self.USER[variant_id].keys()))};")
                    
                    f_out.write(f"CDSS={score}\n")
            
            # Non-coding variants
            for score in sorted(nc_score.keys(), reverse=True):
                for variant_id in sorted(nc_score[score].keys()):
                    if variant_id in self.VCF:
                        f_out.write(f"{self.VCF[variant_id]}\t")
                    elif variant_id in self.DES:
                        fields = self.DES[variant_id].split('\t')
                        if len(fields) >= 5:
                            f_out.write(f"{fields[0]}\t{fields[2]}\t.\t{fields[3]}\t{fields[4]}\t.\t.\t")
                    else:
                        continue
                    
                    f_out.write(f"SAMPLE={sample};")
                    id_short = '\t'.join(variant_id.split('\t')[:2])
                    f_out.write(f"GERP={self.GERP.get(id_short, '.')};")
                    f_out.write("CDS=No;")
                    
                    if variant_id in self.HUB:
                        f_out.write(f"HUB={','.join(sorted(self.HUB[variant_id].keys()))};")
                    if variant_id in self.ANNO:
                        f_out.write(f"NCENC={','.join(sorted(self.ANNO[variant_id].keys()))};")
                    if variant_id in self.HOT:
                        f_out.write(f"HOT={','.join(sorted(self.HOT[variant_id].keys()))};")
                    if variant_id in self.MOTIFBR:
                        f_out.write(f"MOTIFBR={self.MOTIFBR[variant_id]};")
                    if variant_id in self.MOTIFG:
                        f_out.write(f"MOTIFG={self.MOTIFG[variant_id]};")
                    if variant_id in self.SEN:
                        f_out.write("SEN=Yes;")
                    if variant_id in self.USEN:
                        f_out.write("USEN=Yes;")
                    if variant_id in self.CONS:
                        f_out.write("UCONS=Yes;")
                    
                    if variant_id in self.GENE:
                        gene_info_str = ""
                        cang_str = ""
                        for gene in sorted(self.GENE[variant_id].keys()):
                            info = f"{gene}({'&'.join(sorted(self.GENE[variant_id][gene].keys()))})"
                            gene_info_str += info
                            if gene in gene_info:
                                cang_str += f"{gene}{''.join(sorted(gene_info[gene].keys()))},"
                        
                        f_out.write(f"GENE={gene_info_str};")
                        if cang_str:
                            f_out.write(f"CANG={cang_str.rstrip(',')};")
                    
                    if hasattr(self, 'USER') and variant_id in self.USER:
                        f_out.write(f"USER_ANNO={','.join(sorted(self.USER[variant_id].keys()))};")
                    
                    f_out.write(f"NCDS={score}\n")
        
        with open(out, 'w') as f_out:
            if output_format.lower() == 'bed':
                print_bed_format()
            elif output_format.lower() == 'vcf':
                print_vcf_format()
    
    @staticmethod
    def recur(out_detail: str, out_recur: str, out_driver: str, outformat: str,
              cancer_dir: str, cancer_type: str, score_cut: float, anno_dir: str,
              weight_file: str, recur_db_use: str):
        """Recurrent analysis (recurrent in the same position, same gene or same regulatory elements)."""
        # Complex recurrence analysis function
        # This would analyze recurrence across samples
        # Placeholder implementation
        pass  # Full implementation needed

