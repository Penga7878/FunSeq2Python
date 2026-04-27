#!/usr/bin/env python3
"""
FunSeq2 Python Implementation
A flexible framework to annotate and prioritize cancer somatic mutations.
"""

import os
import sys
import re
import argparse
import subprocess
import multiprocessing
import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Import our modules
from Funseq_SNV import Funseq_SNV
from Funseq_Indel import Funseq_Indel

'''
#Reads a text context file
#looks for file_path =...
#looks for key=value lines
returns a dictionary of the key=value pairs
(
 "/home/marcus/data",
 {
   "input": "/home/marcus/data/data.csv",
   "output": "/home/marcus/data/result.txt"
 }
)
'''
def read_config(config_file: str) -> Tuple[str, Dict[str, str]]:
    """Read configuration file and return file_path and variables dictionary."""
    file_path = ""
    variables = {}
    
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if line.startswith('file_path='):
                file_path = line.split('=', 1)[1]
            elif '=' in line:
                match = re.match(r'^(\w+)\s*=\s*(.+)$', line)
                if match:
                    key, value = match.groups()
                    variables[key] = os.path.join(file_path, value)
    
    return file_path, variables

'''
validates user inputs by checking ranges and allowed values for parameters like MAF, genome mode, file formats, and cutoffs
'''
def validate_parameters(maf: float, genome_mode: int, informat: str, 
                       outformat: str, sv_length_cut: str, exp_format: str = ""):
    """Validate input parameters."""
    if not (0 <= maf <= 1):
        sys.exit("Error. Please enter proper MAF...\n")
    
    if genome_mode not in [1, 2]:
        sys.exit("Error. Please specify correct Genome Mode...\n")
    
    if informat.lower() not in ['bed', 'vcf']:
        sys.exit("Error. Input format should be bed or vcf...\n")
    
    if outformat.lower() not in ['bed', 'vcf']:
        sys.exit("Error. Output format should be bed or vcf...\n")
    
    if sv_length_cut != 'inf' and not sv_length_cut.isdigit():
        sys.exit("Error. Indel length cutoff should be integer or 'inf'...\n")
    
    if exp_format and exp_format.lower() not in ['rpkm', 'raw']:
        sys.exit("Error. Expression format should be rpkm or raw...\n")


def check_file(filepath: str, error_msg: str):
    """Check if file exists and is not empty."""
    if not (os.path.isfile(filepath) and os.path.getsize(filepath) > 0):
        sys.exit(f"Error: {error_msg} not found or empty...\n")


def main_analysis(infile: str, output_path: str, tag: str, variables: Dict[str, str],
                  maf: float, informat: str, outformat: str, genome_mode: int,
                  nc_mode: int, sv_length_cut: str, motif_p_value_cut: float,
                  user_anno_dir: str, weight_mode: int, weight_file: str,
                  gene_list: str = "", expression: str = "", class_file: str = "",
                  exp_format: str = ""):
    """Main analysis pipeline for a single input file."""
    
    # Check input file
    if not (os.path.isfile(infile) and os.path.getsize(infile) > 0):
        sys.exit("Error: input file not found or empty...\n")
    
    # Clean up previous outputs
    for pattern in glob.glob(f"{output_path}/{tag}.*"):
        os.remove(pattern)
    
    # Error log file
    err_file = f"{output_path}/{tag}.err"
    sys.stderr = open(err_file, 'w')
    
    # Intermediate files
    out_snv_filter = f"{output_path}/{tag}.out.snv.filter"
    out_nc = f"{output_path}/{tag}.out.nc"
    out_cds = f"{output_path}/{tag}.out.cds"
    out_vat = f"{output_path}/{tag}.out.vat"
    out_motif = f"{output_path}/{tag}.out.motif"
    out_indel = f"{output_path}/{tag}.out.indel"
    
    # Output file
    run_out = f"{output_path}/{tag}.result.{outformat}"
    
    # Gene list filtering
    '''If the user provides a gene list, this code:

        #1. Finds all genomic regions (CDS, promoters, introns, UTRs, enhancers) that belong to those genes.

        #2. Combines them into one BED file.

        #3. Filters the input variants to keep only variants that overlap those regions
    '''
    if gene_list:
        out_gene_sel = f"{output_path}/{tag}.out.gene.sel"
        cds = variables['cds']
        promoter = variables['promoter']
        intron = variables['intron']
        utr = variables['utr']
        enhancer = variables['enhancer']
        
        
        # Extract regions for genes in gene_list
        subprocess.run(f"awk '{{print \"\\t\"$0\"$\"}}' {gene_list} | grep -f - {cds} | cut -f 1,2,3 > {out_gene_sel}", shell=True)
        subprocess.run(f"awk '{{print \"\\t\"$0\"$\"}}' {gene_list} | grep -f - {promoter} | cut -f 1,2,3 >> {out_gene_sel}", shell=True)
        subprocess.run(f"awk '{{print \"\\t\"$0\"$\"}}' {gene_list} | grep -f - {intron} | cut -f 1,2,3 >> {out_gene_sel}", shell=True)
        subprocess.run(f"awk '{{print \"\\t\"$0\"$\"}}' {gene_list} | grep -f - {utr} | cut -f 1,2,3 >> {out_gene_sel}", shell=True)
        subprocess.run(f"awk '{{print \"\\t\"$0\"$\"}}' {gene_list} | grep -f - {enhancer} | cut -f 1,2,3 >> {out_gene_sel}", shell=True)
        
        filtered_infile = f"{output_path}/{tag}.infile"
        subprocess.run(f"intersectBed -u -a {infile} -b {out_gene_sel} > {filtered_infile}", shell=True)
        infile = filtered_infile
        
        if not (os.path.isfile(infile) and os.path.getsize(infile) > 0):
            sys.exit("Error: no variants occurred in / associated with requested genes...\n")
    
    print(f"------------- Running: {tag} starts (0%)---------------\n")
    
    de_data = f"{output_path}/DE.gene.txt"
    
    # Format check: checking that your input variant file (infile) matches the declared format (informat)
    print(f"... Input format check : {informat} ...\n")
    status = Funseq_SNV.format_check(infile, informat, err_file)
    if status == 0:
        print("... Format ok ...\n")
    else:
        sys.exit(1)
    
    # Step 1: 1000 genomes SNV filtering
    print(f"... Start filtering SNVs with minor allele frequency = {maf} ...\n")
    
    data = Funseq_SNV()
    data.snv_filter(infile, informat, variables['tgp_snp'], maf, out_snv_filter, 
                   out_indel, sv_length_cut)
    #deletes temporary file if exists
    if os.path.exists(f"{output_path}/{tag}.infile"):
        os.remove(f"{output_path}/{tag}.infile")

    #Checks if every SNV matched common population variants
    if os.path.getsize(out_snv_filter) == 0:
        print(f"Warning: sample {tag} - no SNVs left after filtering against natural variations ...\n")
        sys.stderr.write("Warning: no SNVs left after filtering against natural variations ...\n")
    else:
        #Otherwise: add GERP conservation scores to the SNVs
        data.gerp_score(variables['gerp_file'], out_snv_filter)
        

        '''
        Take the filtered SNVs (out_snv_filter) and:

        #1. split them into coding (CDS-overlapping) vs non-coding variants

        #2. annotate the non-coding ones with ENCODE-style features

        #3. apply any user-provided annotations to all filtered SNVs
        '''
        # Step 2: Get Annotated (Coding & Non-coding) SNVs
        print("... Separate into non-coding annotation & Coding (30%)...\n")
        
        subprocess.run(f"intersectBed -u -a {out_snv_filter} -b {variables['cds']} > {out_cds}", shell=True)
        if os.path.getsize(out_cds) > 0:
            subprocess.run(f"intersectBed -v -a {out_snv_filter} -b {out_cds} > {out_nc}", shell=True)
        else:
            subprocess.run(f"cp {out_snv_filter} {out_nc}", shell=True)
        
        if os.path.getsize(out_nc) > 0:
            data.read_encode(out_nc, variables['encode_annotation'])
        
        data.user_annotation(out_snv_filter, user_anno_dir)
        
        # Step 3: Non-coding
        print("... Start Non-coding SNVs analysis (50%)...\n")
        
        data.conserved(out_snv_filter, variables['conserved'])
        if os.path.getsize(out_nc) > 0:
            data.sensitive(out_nc, variables['sensitive'])
            data.hot_region(out_nc, variables['hot_file'])
            data.motif_break(out_nc, variables['ancestral_file'], variables['bound_motif'],
                           genome_mode, variables['motif_pfm'])
            data.gene_link(out_nc, variables['promoter'], variables['enhancer'],
                          variables['intron'], variables['utr'], variables['network_dir'])
            if variables.get('reference_file'):
                data.motif_gain(variables['motif_pfm'], variables['reference_file'],
                              variables['score_file'], motif_p_value_cut, out_motif)
        
        # Step 4: Coding
        print("... Start Coding SNVs analysis (80%)...\n")
        
        if os.path.getsize(out_cds) > 0:
            status = Funseq_SNV.coding(out_cds, variables['coding_interval'],
                                     variables['coding_fasta'], variables['network_dir'],
                                     variables['selection'], out_vat, output_path,
                                     nc_mode, variables['cds'])
            if status == 1:
                print("Warning: coding variants not analyzed by VAT ... \n")
                sys.stderr.write("Warning: coding variants not analyzed by VAT ... \n")
        
        # Step 5: Scoring SNVs
        data.intergrate(outformat, tag, out_nc, out_vat, variables['gene_info_dir'],
                       variables['reg_net'], run_out, de_data, weight_mode, weight_file)

    # Indel analysis
    # `out_indel` is only created when there are indels; guard against missing file.
    if os.path.exists(out_indel) and os.path.getsize(out_indel) > 0:
        print("... Start Indels analysis ...\n")
        
        indel = Funseq_Indel()
        run_indel_out = f"{output_path}/{tag}.indel.{outformat}"
        out_indel_motif = f"{output_path}/{tag}.out.indel.motif"
        
        indel.gerp_score(variables['gerp_file'], out_indel)
        indel.annotations(out_indel, variables['conserved'], variables['hot_file'],
                         variables['sensitive'], variables['encode_annotation'],
                         user_anno_dir, variables['bound_motif'])
        indel.gene_link(out_indel, variables['promoter'], variables['enhancer'],
                       variables['intron'], variables['utr'], variables['network_dir'],
                       variables['cds'], variables['selection'])
        indel.coding(out_indel, variables['coding_interval'], variables['coding_fasta'], nc_mode)
        if variables.get('reference_file'):
            indel.motif_gain(variables['motif_pfm'], variables['reference_file'],
                            variables['score_file'], motif_p_value_cut, out_indel_motif)
        indel.intergrate(outformat, tag, variables['gene_info_dir'], variables['reg_net'],
                        run_indel_out, de_data)
    
    print(f"------------- {tag} Done (100%)---------------\n")
    
    # Cleanup intermediate files (optional)
    choice = 0  # 1 will delete all intermediate results
    if choice == 1:
        for pattern in glob.glob(f"{output_path}/{tag}.out*"):
            os.remove(pattern)
    
    sys.stderr.close()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='FunSeq2: Annotate and prioritize cancer somatic mutations')
    
    parser.add_argument('infile', help='Input variants file')
    parser.add_argument('maf', type=float, help='Minor allele frequency threshold')
    parser.add_argument('genome_mode', type=int, choices=[1, 2], 
                       help='1: somatic; 2: germline | personal genome')
    parser.add_argument('informat', choices=['bed', 'vcf'], help='Input format [bed | vcf]')
    parser.add_argument('outformat', choices=['bed', 'vcf'], help='Output format [bed | vcf]')
    parser.add_argument('nc_mode', type=int, help='Only do non-coding if nc_mode=1')
    parser.add_argument('output_path', help='Output path')
    parser.add_argument('num_per_run', type=int, help='Number of genome per run')
    parser.add_argument('cancer_type', help='Cancer type retrieved')
    parser.add_argument('score_cut', type=float, help='Non-coding candidate score cut')
    parser.add_argument('weight_mode', type=int, choices=[0, 1],
                       help='0: unweighted; 1: weighted scoring scheme')
    parser.add_argument('user_anno_dir', help='Directory for user-specific annotations')
    parser.add_argument('recur_db_use', help='Recurrence database use flag')
    parser.add_argument('sv_length_cut', help='SV length cut')
    parser.add_argument('gene_list', nargs='?', default='', help='Gene list file (optional)')
    parser.add_argument('--expression', default='', help='Expression file (optional)')
    parser.add_argument('--class', dest='class_file', default='', help='Sample class file (optional)')
    parser.add_argument('--exp_format', default='', help='Expression format: rpkm/raw (optional)')
    parser.add_argument('--config', default='config.txt', help='Configuration file path')
    
    args = parser.parse_args()
    
    # Handle expression/class/exp_format if provided
    if args.expression and args.class_file:
        exp_format = args.exp_format if args.exp_format else 'rpkm'
    else:
        exp_format = ''
    
    # Validate parameters
    validate_parameters(args.maf, args.genome_mode, args.informat, args.outformat,
                       args.sv_length_cut, exp_format)
    
    # Create output directory
    os.makedirs(args.output_path, exist_ok=True)
    
    # Read configuration
    config_file = os.environ.get('FunSeqConfig', args.config)
    print(f"FunSeq is starting!\n############################\nConfig file:{config_file}\n")
    
    file_path, variables = read_config(config_file)
    
    # Set up GENCODE paths
    gencode_dir = variables.get('gencode', '')
    if gencode_dir:
        # Find gencode version
        result = subprocess.run(f"ls {gencode_dir}/*.promoter.bed", shell=True, 
                              capture_output=True, text=True)
        if result.returncode == 0 and result.stdout:
            gencode_file = result.stdout.strip().split('\n')[-1]
            gencode_v = Path(gencode_file).stem.split('.')[1]
            variables['cds'] = f"{gencode_dir}/gencode.{gencode_v}.cds.bed"
            variables['promoter'] = f"{gencode_dir}/gencode.{gencode_v}.promoter.bed"
            variables['intron'] = f"{gencode_dir}/gencode.{gencode_v}.intron.bed"
            variables['utr'] = f"{gencode_dir}/gencode.{gencode_v}.utr.bed"
            variables['coding_interval'] = f"{gencode_dir}/gencode.{gencode_v}.cds.interval"
            variables['coding_fasta'] = f"{gencode_dir}/gencode.{gencode_v}.cds.fa"
    
    # Check required files
    required_files = {
        'tgp_snp': variables.get('tgp_snp'),
        'encode_annotation': variables.get('encode_annotation'),
        'cds': variables.get('cds'),
        'bound_motif': variables.get('bound_motif'),
        'promoter': variables.get('promoter'),
        'enhancer': variables.get('enhancer'),
        'motif_pfm': variables.get('motif_pfm'),
        'score_file': variables.get('score_file'),
        'weight_file': variables.get('weight_file'),
    }
    
    for key, filepath in required_files.items():
        if filepath:
            check_file(filepath, key)
    
    # Reference file is only needed for motif_gain, check it separately
    reference_file = variables.get('reference_file')
    if reference_file and not (os.path.isfile(reference_file) and os.path.getsize(reference_file) > 0):
        print(f"Warning: reference_file ({reference_file}) not found. Motif gain analysis will be skipped.\n")
        variables['reference_file'] = None
    
    if args.genome_mode == 2:
        check_file(variables.get('ancestral_file', ''), 'ancestral_file')
    
    if args.nc_mode == 0:
        check_file(variables.get('coding_interval', ''), 'coding_interval')
        check_file(variables.get('coding_fasta', ''), 'coding_fasta')
    
    # Main output files
    file_detail = f"{args.output_path}/Output.{args.outformat}"
    file_recur = f"{args.output_path}/Recur.Summary"
    file_driver = f"{args.output_path}/Candidates.Summary"
    file_err = f"{args.output_path}/Error.log"
    file_indel = f"{args.output_path}/Output.indel.{args.outformat}"
    
    # Remove old output files
    for f in [file_detail, file_recur, file_driver, file_err, file_indel]:
        if os.path.exists(f):
            os.remove(f)
    
    # Expression analysis
    if args.expression and args.class_file:
        print("Differential Gene Expression Analysis ...\n")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script = os.path.join(os.path.dirname(script_dir), 'differential_gene_expression.r')
        subprocess.run(['Rscript', r_script, args.expression, args.class_file,
                       exp_format, args.output_path, file_err])
    
    # Process input files
    input_files = args.infile.split(',')
    motif_p_value_cut = 4e-8
    
    # Simplified processing (single file case)
    # TODO: Add parallel processing support
    if len(input_files) == 1:
        tag = Path(input_files[0]).stem.split('.')[0]
        main_analysis(input_files[0], args.output_path, tag, variables,
                     args.maf, args.informat, args.outformat, args.genome_mode,
                     args.nc_mode, args.sv_length_cut, motif_p_value_cut,
                     args.user_anno_dir, args.weight_mode, variables['weight_file'],
                     args.gene_list, args.expression, args.class_file, exp_format)
        
        # Move/concatenate output files
        tmp_err = f"{args.output_path}/{tag}.err"
        if os.path.exists(tmp_err):
            os.rename(tmp_err, file_err)
        
        tmp_in = f"{args.output_path}/{tag}.result.{args.outformat}"
        if os.path.exists(tmp_in):
            os.rename(tmp_in, file_detail)
        
        tmp_indel = f"{args.output_path}/{tag}.indel.{args.outformat}"
        if os.path.exists(tmp_indel):
            os.rename(tmp_indel, file_indel)
    
    # Recurrence analysis
    if os.path.exists(file_detail) and os.path.getsize(file_detail) > 0:
        Funseq_SNV.recur(file_detail, file_recur, file_driver, args.outformat,
                        variables['cancer_dir'], args.cancer_type, args.score_cut,
                        args.user_anno_dir, variables['weight_file'], args.recur_db_use)
    
    if os.path.exists(file_indel) and os.path.getsize(file_indel) > 0:
        Funseq_Indel.combine(file_indel, args.outformat)


if __name__ == '__main__':
    main()

