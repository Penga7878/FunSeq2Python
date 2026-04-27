"""
Funseq_Indel Python Module
Handles Indel (Insertion/Deletion) analysis
"""

import os
import re
import subprocess
from typing import Dict, List, Optional


class Funseq_Indel:
    """Class for Indel analysis and annotation."""
    
    def __init__(self):
        """Initialize the object with empty dictionaries for data storage."""
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
        self.MOTIFG = {}  # Motif gain
        self.VAT = {}  # VAT annotations
        self.SELECTION = {}  # Selection data
    
    def gerp_score(self, gerp_file: str, infile: str):
        """Obtain GERP++ scores for indels."""
        # Initialize GERP scores
        with open(infile, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 5:
                    variant_id = '\t'.join(fields[:5])
                    self.GERP[variant_id] = "."
        
        if os.path.isfile(gerp_file) and os.path.getsize(gerp_file) > 0:
            cmd = f"awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print $1,$2,$3,$1\":\"$2\":\"$3\":\"$4\":\"$5}}' {infile} | sort -k 1,1 -k 2,2n | uniq | bigWigAverageOverBed {gerp_file} stdin stdout"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = re.split(r'[:|\s]+', line)
                    if len(parts) >= 6:
                        variant_id = '\t'.join(parts[:5])
                        score = float(parts[-1])
                        if score != 0:
                            self.GERP[variant_id] = score
    
    def annotations(self, input_file: str, conserved: str, hot: str, sensitive: str,
                    encode: str, anno_dir: str, bound_motif: str):
        """Add annotations for indels."""
        # Conserved regions
        cmd = f"intersectBed -u -a {input_file} -b {conserved}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = '\t'.join(fields[:5])
                    self.CONS[variant_id] = 1
        
        # HOT regions
        cmd = f"intersectBed -a {input_file} -b {hot} -wo"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 9:
                    variant_id = '\t'.join(fields[:5])
                    if variant_id not in self.HOT:
                        self.HOT[variant_id] = {}
                    self.HOT[variant_id][fields[8]] = 1
        
        # Sensitive regions
        cmd = f"intersectBed -u -a {input_file} -b {sensitive}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = '\t'.join(fields[:5])
                    self.SEN[variant_id] = 1
        
        # Ultra-sensitive regions
        cmd = f"grep 'Ultra' {sensitive} | intersectBed -u -a {input_file} -b stdin"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split()
                if len(fields) >= 5:
                    variant_id = '\t'.join(fields[:5])
                    self.USEN[variant_id] = 1
        
        # ENCODE annotations
        cmd = f"intersectBed -a {input_file} -b {encode} -wo"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split('\t')
                if len(fields) >= 9:
                    variant_id = '\t'.join(fields[:5])
                    interval = f"{fields[5]}:{fields[6]}-{fields[7]}"
                    anno_parts = fields[8].split('.')
                    
                    if variant_id not in self.ANNO:
                        self.ANNO[variant_id] = {}
                    
                    if anno_parts and anno_parts[-1] == 'u':
                        tag = f"{anno_parts[0]}({'.'.join(anno_parts[1:-1])})"
                    else:
                        tag = f"{anno_parts[0]}({'.'.join(anno_parts[1:])}|{interval})"
                    self.ANNO[variant_id][tag] = 1
        
        # User annotations
        if os.path.isdir(anno_dir):
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
                                variant_id = '\t'.join(fields[:5])
                                interval = f"{fields[5]}:{fields[6]}-{fields[7]}"
                                
                                if variant_id not in self.USER:
                                    self.USER[variant_id] = {}
                                
                                if len(fields) > 9:
                                    tag = f"{cate}({fields[8]}|{interval})"
                                else:
                                    tag = f"{cate}({interval})"
                                self.USER[variant_id][tag] = 1
        
        # Bound motif annotations
        cmd = f"intersectBed -a {bound_motif} -b {input_file} -wo | sort -k 1,1 -k 2,2n | uniq"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split('\t')
                if len(fields) >= 12:
                    variant_id = '\t'.join(fields[7:12])  # Adjusted indices
                    interval = f"{fields[0]}:{fields[1]}-{fields[2]}"
                    
                    if variant_id not in self.ANNO:
                        self.ANNO[variant_id] = {}
                    
                    tag = f"TFM({fields[11]}|{fields[8]}|{interval})"
                    self.ANNO[variant_id][tag] = 1
                    
                    # Motif breaking
                    motif_key = f"{fields[11]}#{fields[8]}#{fields[1]}#{fields[2]}#{fields[10]}"
                    if variant_id in self.MOTIFBR:
                        self.MOTIFBR[variant_id] = f"{self.MOTIFBR[variant_id]},{motif_key}"
                    else:
                        self.MOTIFBR[variant_id] = motif_key
    
    def gene_link(self, input_file: str, promoter: str, distal: str, intron: str,
                  utr: str, network: str, cds: str, selection: str):
        """Link indels with genes (promoters & regulatory elements)."""
        # Load network data
        network_data = {}
        net_degree = {}
        if os.path.isdir(network):
            for filename in os.listdir(network):
                if not filename.startswith('.'):
                    filepath = os.path.join(network, filename)
                    net_name = filename.split('.')[0]
                    with open(filepath, 'r') as f:
                        for line in f:
                            if 'GENE_NAME' not in line:
                                fields = line.split()
                                if len(fields) >= 2:
                                    gene = fields[0]
                                    degree = float(fields[1])
                                    if gene not in network_data:
                                        network_data[gene] = {}
                                    network_data[gene][net_name] = degree
                                    
                                    if net_name not in net_degree:
                                        net_degree[net_name] = []
                                    net_degree[net_name].append(degree)
        
        # Load selection data
        selection_genes = {}
        if os.path.isfile(selection):
            with open(selection, 'r') as f:
                for line in f:
                    if 'GENE_NAME' not in line:
                        gene = line.split()[0]
                        selection_genes[gene] = 1
        
        # Intersect with different gene regions
        for cmd, tag in [
            (f"intersectBed -a {input_file} -b {intron} -wo | sort -k 1,1 -k 2,2n | uniq", "Intron"),
            (f"intersectBed -a {input_file} -b {utr} -wo | sort -k 1,1 -k 2,2n | uniq", "UTR"),
            (f"intersectBed -a {input_file} -b {promoter} -wo | sort -k 1,1 -k 2,2n | uniq", "Promoter"),
            (f"intersectBed -a {input_file} -b {distal} -wo | sort -k 1,1 -k 2,2n | uniq", "Distal"),
            (f"intersectBed -a {input_file} -b {cds} -wo | sort -k 1,1 -k 2,2n | uniq", "Coding"),
        ]:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in result.stdout.strip().split('\n'):
                if line:
                    fields = line.split()
                    if len(fields) >= 9:
                        variant_id = '\t'.join(fields[:5])
                        gene = fields[8]
                        
                        if variant_id not in self.GENE:
                            self.GENE[variant_id] = {}
                        if gene not in self.GENE[variant_id]:
                            self.GENE[variant_id][gene] = {}
                        self.GENE[variant_id][gene][tag] = 1
                        
                        # Selection
                        if gene in selection_genes:
                            if variant_id not in self.SELECTION:
                                self.SELECTION[variant_id] = {}
                            self.SELECTION[variant_id][gene] = 1
                        
                        # Network hubs
                        if gene in network_data:
                            prob_genes = []
                            for net_name in sorted(network_data[gene].keys()):
                                degree = network_data[gene][net_name]
                                if net_name in net_degree:
                                    greater = sum(1 for d in net_degree[net_name] if d > degree)
                                    prob = 1 - (greater / len(net_degree[net_name]))
                                    prob_genes.append(f"{net_name}({prob:.3f})")
                            
                            hub_key = f"{gene}:{''.join(prob_genes)}"
                            if variant_id not in self.HUB:
                                self.HUB[variant_id] = {}
                            self.HUB[variant_id][hub_key] = 1
    
    def coding(self, snp_input: str, file_interval: str, file_fasta: str, nc_mode: int):
        """Coding pipe with VAT (variant annotation tool) analysis for indels."""
        if nc_mode == 0:
            cmd = f"awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}}{{print \"##fileformat=VCFv4.0\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\"}}{{print $1,$3,$2,$4,$5,\".\",\"PASS\",\".\"}}' {snp_input} | indelMapper {file_interval} {file_fasta}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            for line in result.stdout.strip().split('\n'):
                if line and not line.startswith('#'):
                    fields = line.split()
                    if len(fields) >= 8:
                        variant_id = f"{fields[0]}\t{fields[2]}\t{fields[1]}\t{fields[3]}\t{fields[4]}"
                        self.VAT[variant_id] = fields[7]
    
    def motif_gain(self, pfm_file: str, reference_file: str, score_file: str,
                   p_cut: float, out_tmp: str):
        """Check whether indels create novel motifs."""
        # Complex function for motif gain analysis
        # Requires sequence extraction, PWM calculation, and scoring
        # This is a placeholder - full implementation would be similar to SNV version
        # but adapted for indels
        pass  # Full implementation needed
    
    def intergrate(self, output_format: str, sample: str, gene_dir: str, reg_net: str,
                   out: str, de_data: str):
        """Integrate all outputs for indels."""
        # Read gene information
        gene_info = {}
        if os.path.isdir(gene_dir):
            for filename in os.listdir(gene_dir):
                if not filename.startswith('.'):
                    filepath = os.path.join(gene_dir, filename)
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
        
        # Output
        with open(out, 'w') as f_out:
            if output_format.lower() == 'bed':
                for variant_id in sorted(self.GERP.keys()):
                    fields = variant_id.split('\t')
                    if len(fields) >= 5:
                        f_out.write(f"{variant_id}\t{sample}\t")
                        f_out.write(f"{self.GERP.get(variant_id, '.')};")
                        f_out.write(f"{self.VAT.get(variant_id, '.')};")
                        
                        # HUB
                        if variant_id in self.HUB:
                            f_out.write(f"{','.join(sorted(self.HUB[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        # SELECTION
                        if variant_id in self.SELECTION:
                            f_out.write(f"{','.join(sorted(self.SELECTION[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        # ANNO
                        if variant_id in self.ANNO:
                            f_out.write(f"{','.join(sorted(self.ANNO[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        # HOT
                        if variant_id in self.HOT:
                            f_out.write(f"{','.join(sorted(self.HOT[variant_id].keys()))};")
                        else:
                            f_out.write(".;")
                        
                        # MOTIF
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
                        
                        # SEN
                        f_out.write("Yes;" if variant_id in self.SEN else ".;")
                        # USEN
                        f_out.write("Yes;" if variant_id in self.USEN else ".;")
                        # CONS
                        f_out.write("Yes;" if variant_id in self.CONS else ".;")
                        
                        # GENE
                        if variant_id in self.GENE:
                            gene_list = []
                            for gene in sorted(self.GENE[variant_id].keys()):
                                info = f"{gene}({'&'.join(sorted(self.GENE[variant_id][gene].keys()))})"
                                if gene in gene_info:
                                    info += ''.join(sorted(gene_info[gene].keys()))
                                gene_list.append(info)
                            f_out.write(f"{','.join(gene_list)};")
                        else:
                            f_out.write(".;")
                        
                        # USER
                        if variant_id in self.USER:
                            f_out.write(f"{','.join(sorted(self.USER[variant_id].keys()))}\n")
                        else:
                            f_out.write(".\n")
            
            elif output_format.lower() == 'vcf':
                for variant_id in sorted(self.GERP.keys()):
                    fields = variant_id.split('\\t')
                    if len(fields) >= 5:
                        f_out.write(f"{fields[0]}\t{int(fields[1])+1}\t.\t{fields[3]}\t{fields[4]}\t.\t.\t")
                        f_out.write(f"SAMPLE={sample};")
                        f_out.write(f"GERP={self.GERP.get(variant_id, '.')};")
                        
                        if variant_id in self.VAT:
                            f_out.write(f"{self.VAT[variant_id]};")
                        
                        if variant_id in self.HUB:
                            f_out.write(f"HUB={','.join(sorted(self.HUB[variant_id].keys()))};")
                        
                        if variant_id in self.SELECTION:
                            f_out.write(f"GNEG={','.join(sorted(self.SELECTION[variant_id].keys()))};")
                        
                        if variant_id in self.ANNO:
                            f_out.write(f"NCENC={','.join(sorted(self.ANNO[variant_id].keys()))};")
                        
                        if variant_id in self.HOT:
                            f_out.write(f"HOT={','.join(sorted(self.HOT[variant_id].keys()))};")
                        
                        if variant_id in self.MOTIFBR:
                            f_out.write(f"MOTIFBR={self.MOTIFBR[variant_id]};")
                        
                        if variant_id in self.MOTIFG:
                            f_out.write(f"MOTIFG={self.MOTIFG[variant_id]};")
                        
                        f_out.write("SEN=Yes;" if variant_id in self.SEN else "")
                        f_out.write("USEN=Yes;" if variant_id in self.USEN else "")
                        f_out.write("UCONS=Yes;" if variant_id in self.CONS else "")
                        
                        if variant_id in self.GENE:
                            gene_list = []
                            gene_info_str = ""
                            for gene in sorted(self.GENE[variant_id].keys()):
                                info = f"{gene}({'&'.join(sorted(self.GENE[variant_id][gene].keys()))})"
                                gene_list.append(info)
                                if gene in gene_info:
                                    gene_info_str += f"{gene}{''.join(sorted(gene_info[gene].keys()))},"
                            
                            f_out.write(f"GENE={','.join(gene_list)};")
                            if gene_info_str:
                                f_out.write(f"CANG={gene_info_str.rstrip(',')};")
                        
                        if variant_id in self.USER:
                            f_out.write(f"USER_ANNO={','.join(sorted(self.USER[variant_id].keys()))}")
                        
                        f_out.write("\n")
    
    @staticmethod
    def combine(file: str, format_type: str):
        """Add header to output file."""
        if format_type.lower() == 'bed':
            header = "\t".join(["#chr", "start", "end", "ref", "alt", "sample",
                                "gerp;variant.annotation.cds;network.hub;gene.under.negative.selection;ENCODE.annotated;hot.region;motif.analysis;sensitive;ultra.sensitive;ultra.conserved;target.gene[known_cancer_gene/TF_regulating_known_cancer_gene,differential_expressed_in_cancer,actionable_gene];user.annotations"])
            subprocess.run(f"sed -i '1i {header}' {file}", shell=True)
        else:
            header = "\n".join([
                "##fileformat=VCFv4.0",
                '##INFO=<ID=SAMPLE,Number=.,Type=String,Description="Sample id">',
                '##INFO=<ID=VA,Number=.,Type=String,Description="Coding Variant Annotation">',
                '##INFO=<ID=HUB,Number=.,Type=String,Description="Network Hubs...">',
                '##INFO=<ID=GNEG,Number=.,Type=String,Description="Gene Under Negative Selection">',
                '##INFO=<ID=GERP,Number=.,Type=String,Description="Gerp Score">',
                '##INFO=<ID=NCENC,Number=.,Type=String,Description="NonCoding ENCODE Annotation">',
                '##INFO=<ID=HOT,Number=.,Type=String,Description="Highly Occupied Target Region">',
                '##INFO=<ID=MOTIFBR,Number=.,Type=String,Description="Motif Breaking">',
                '##INFO=<ID=MOTIFG,Number=.,Type=String,Description="Motif Gain">',
                '##INFO=<ID=SEN,Number=.,Type=String,Description="In Sensitive Region">',
                '##INFO=<ID=USEN,Number=.,Type=String,Description="In Ultra-Sensitive Region">',
                '##INFO=<ID=UCONS,Number=.,Type=String,Description="In Ultra-Conserved Region">',
                '##INFO=<ID=GENE,Number=.,Type=String,Description="Target Gene...">',
                '##INFO=<ID=CANG,Number=.,Type=String,Description="Prior Gene Information...">',
                '##INFO=<ID=USER_ANNO,Number=.,Type=String,Description="Annotations from user-input">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            ])
            subprocess.run(f"sed -i '1i {header}' {file}", shell=True)

