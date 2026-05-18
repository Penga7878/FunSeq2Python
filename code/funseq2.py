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
import shlex
import multiprocessing
import glob
from pathlib import Path
from typing import Dict, List, Tuple

from Funseq_SNV import Funseq_SNV
from Funseq_Indel import Funseq_Indel


def read_config(config_file: str) -> Tuple[str, Dict[str, str]]:
    """Read configuration file and return file_path and variables dictionary."""
    file_path = ""
    variables = {}
    lines = []

    config_dir = os.path.dirname(os.path.abspath(config_file))

    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            lines.append(line)

    # First pass: extract file_path so later keys can be joined correctly
    for line in lines:
        if line.startswith('file_path='):
            file_path = line.split('=', 1)[1].strip()
            break

    # Resolve relative file_path against config file's directory
    if file_path and not os.path.isabs(file_path):
        file_path = os.path.join(config_dir, file_path)

    # Second pass: build variables dict
    for line in lines:
        if line.startswith('file_path='):
            continue
        match = re.match(r'^(\w+)\s*=\s*(.+)$', line)
        if match:
            key, value = match.groups()
            variables[key] = os.path.join(file_path, value.strip())

    return file_path, variables


def validate_parameters(maf: float, genome_mode: int, informat: str,
                        outformat: str, sv_length_cut: str, exp_format: str = ""):
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
    if not (os.path.isfile(filepath) and os.path.getsize(filepath) > 0):
        sys.exit(f"Error: {error_msg} not found or empty...\n")


def file_nonempty(filepath: str) -> bool:
    return os.path.isfile(filepath) and os.path.getsize(filepath) > 0


def main_analysis(infile: str, output_path: str, tag: str, variables: Dict[str, str],
                  maf: float, informat: str, outformat: str, genome_mode: int,
                  nc_mode: int, sv_length_cut: str, motif_p_value_cut: float,
                  user_anno_dir: str, weight_mode: int, weight_file: str,
                  gene_list: str = "", expression: str = "", class_file: str = "",
                  exp_format: str = ""):
    """Main analysis pipeline for a single input file."""

    if not file_nonempty(infile):
        sys.exit("Error: input file not found or empty...\n")

    for pattern in glob.glob(f"{output_path}/{tag}.*"):
        os.remove(pattern)

    err_file = f"{output_path}/{tag}.err"
    out_snv_filter = f"{output_path}/{tag}.out.snv.filter"
    out_nc = f"{output_path}/{tag}.out.nc"
    out_cds = f"{output_path}/{tag}.out.cds"
    out_vat = f"{output_path}/{tag}.out.vat"
    out_motif = f"{output_path}/{tag}.out.motif"
    out_indel = f"{output_path}/{tag}.out.indel"
    run_out = f"{output_path}/{tag}.result.{outformat}"

    if gene_list:
        out_gene_sel = f"{output_path}/{tag}.out.gene.sel"
        q_gene = shlex.quote(gene_list)
        q_out = shlex.quote(out_gene_sel)
        for i, src_key in enumerate(['cds', 'promoter', 'intron', 'utr', 'enhancer']):
            redirect = '>' if i == 0 else '>>'
            q_src = shlex.quote(variables[src_key])
            subprocess.run(
                f"awk '{{print \"\\t\"$0\"$\"}}' {q_gene} | grep -f - {q_src} | cut -f 1,2,3 {redirect} {q_out}",
                shell=True
            )

        filtered_infile = f"{output_path}/{tag}.infile"
        subprocess.run(
            f"intersectBed -u -a {shlex.quote(infile)} -b {q_out} > {shlex.quote(filtered_infile)}",
            shell=True
        )
        infile = filtered_infile
        if not file_nonempty(infile):
            sys.exit("Error: no variants occurred in / associated with requested genes...\n")

    print(f"------------- Running: {tag} starts (0%)---------------\n")

    de_data = f"{output_path}/DE.gene.txt"

    original_stderr = sys.stderr
    err_fh = open(err_file, 'w')
    sys.stderr = err_fh

    try:
        print(f"... Input format check : {informat} ...\n")
        status = Funseq_SNV.format_check(infile, informat, err_file)
        if status == 0:
            print("... Format ok ...\n")
        else:
            sys.exit(1)

        print(f"... Start filtering SNVs with minor allele frequency = {maf} ...\n")

        data = Funseq_SNV()
        data.snv_filter(infile, informat, variables['tgp_snp'], maf, out_snv_filter,
                        out_indel, sv_length_cut)

        if os.path.exists(f"{output_path}/{tag}.infile"):
            os.remove(f"{output_path}/{tag}.infile")

        if not file_nonempty(out_snv_filter):
            print(f"Warning: sample {tag} - no SNVs left after filtering against natural variations ...\n")
            err_fh.write("Warning: no SNVs left after filtering against natural variations ...\n")
        else:
            data.gerp_score(variables['gerp_file'], out_snv_filter)

            print("... Separate into non-coding annotation & Coding (30%)...\n")

            subprocess.run(
                f"intersectBed -u -a {shlex.quote(out_snv_filter)} -b {shlex.quote(variables['cds'])} > {shlex.quote(out_cds)}",
                shell=True
            )
            if file_nonempty(out_cds):
                subprocess.run(
                    f"intersectBed -v -a {shlex.quote(out_snv_filter)} -b {shlex.quote(out_cds)} > {shlex.quote(out_nc)}",
                    shell=True
                )
            else:
                subprocess.run(f"cp {shlex.quote(out_snv_filter)} {shlex.quote(out_nc)}", shell=True)

            if file_nonempty(out_nc):
                data.read_encode(out_nc, variables['encode_annotation'])

            data.user_annotation(out_snv_filter, user_anno_dir)

            print("... Start Non-coding SNVs analysis (50%)...\n")

            data.conserved(out_snv_filter, variables['conserved'])
            if file_nonempty(out_nc):
                data.sensitive(out_nc, variables['sensitive'])
                data.hot_region(out_nc, variables['hot_file'])
                data.motif_break(out_nc, variables['ancestral_file'], variables['bound_motif'],
                                 genome_mode, variables['motif_pfm'])
                data.gene_link(out_nc, variables['promoter'], variables['enhancer'],
                               variables['intron'], variables['utr'], variables['network_dir'])
                if variables.get('reference_file'):
                    data.motif_gain(variables['motif_pfm'], variables['reference_file'],
                                    variables['score_file'], motif_p_value_cut, out_motif)

            print("... Start Coding SNVs analysis (80%)...\n")

            if file_nonempty(out_cds):
                status = Funseq_SNV.coding(out_cds, variables['coding_interval'],
                                           variables['coding_fasta'], variables['network_dir'],
                                           variables['selection'], out_vat, output_path,
                                           nc_mode, variables['cds'])
                if status == 1:
                    print("Warning: coding variants not analyzed by VAT ... \n")
                    err_fh.write("Warning: coding variants not analyzed by VAT ... \n")

            data.intergrate(outformat, tag, out_nc, out_vat, variables['gene_info_dir'],
                            variables['reg_net'], run_out, de_data, weight_mode, weight_file)

    finally:
        sys.stderr = original_stderr
        err_fh.close()

    if file_nonempty(out_indel):
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


def _process_sample(args_tuple):
    """Worker for multiprocessing.Pool — unpacks into main_analysis."""
    return main_analysis(*args_tuple)


def _collect_results(tags: List[str], output_path: str, outformat: str,
                     file_detail: str, file_err: str, file_indel: str):
    """Append per-tag output files into the combined output files, then remove them."""
    for tag in tags:
        tmp_err = f"{output_path}/{tag}.err"
        tmp_in = f"{output_path}/{tag}.result.{outformat}"
        tmp_indel = f"{output_path}/{tag}.indel.{outformat}"

        if os.path.exists(tmp_err):
            with open(tmp_err) as src, open(file_err, 'a') as dst:
                for line in src:
                    dst.write(f"{tag}\t{line}")
            os.remove(tmp_err)

        if os.path.exists(tmp_in):
            with open(tmp_in) as src, open(file_detail, 'a') as dst:
                dst.write(src.read())
            os.remove(tmp_in)

        if os.path.exists(tmp_indel):
            with open(tmp_indel) as src, open(file_indel, 'a') as dst:
                dst.write(src.read())
            os.remove(tmp_indel)


def process_files(input_files: List[str], output_path: str, variables: Dict[str, str],
                  maf: float, informat: str, outformat: str, genome_mode: int,
                  nc_mode: int, sv_length_cut: str, motif_p_value_cut: float,
                  user_anno_dir: str, weight_mode: int, weight_file: str,
                  gene_list: str, expression: str, class_file: str, exp_format: str,
                  num_per_run: int, file_detail: str, file_err: str, file_indel: str):
    """Dispatch single-file, multi-sample, or multi-file analysis."""

    # args shared across all main_analysis calls (everything after tag)
    common = (variables, maf, informat, outformat, genome_mode, nc_mode,
              sv_length_cut, motif_p_value_cut, user_anno_dir, weight_mode,
              weight_file, gene_list, expression, class_file, exp_format)

    if len(input_files) == 1 and informat.lower() == 'bed':
        result = subprocess.run(
            f"cut -f 6 {shlex.quote(input_files[0])} | sort | uniq",
            shell=True, capture_output=True, text=True
        )
        samples = [s for s in result.stdout.strip().split('\n') if s]
    else:
        samples = []

    if len(input_files) == 1 and len(samples) <= 1:
        # Simple single-sample case
        tag = Path(input_files[0]).stem.split('.')[0] or Path(input_files[0]).stem
        main_analysis(input_files[0], output_path, tag, *common)

        tmp_err = f"{output_path}/{tag}.err"
        if os.path.exists(tmp_err):
            os.rename(tmp_err, file_err)

        tmp_in = f"{output_path}/{tag}.result.{outformat}"
        if os.path.exists(tmp_in):
            os.rename(tmp_in, file_detail)

        tmp_indel = f"{output_path}/{tag}.indel.{outformat}"
        if os.path.exists(tmp_indel):
            os.rename(tmp_indel, file_indel)

    elif len(input_files) == 1 and len(samples) > 1:
        # Single BED with multiple samples identified by column 6
        tags = []
        tasks = []
        for sample_id in samples:
            tag = sample_id.split('.')[0]
            tags.append(tag)
            sample_file = f"{output_path}/{sample_id}"
            subprocess.run(
                f"awk '{{FS=\"\\t\";OFS=\"\\t\"}} $6 == {shlex.quote(sample_id)}' "
                f"{shlex.quote(input_files[0])} > {shlex.quote(sample_file)}",
                shell=True
            )
            tasks.append((sample_file, output_path, tag) + common)

        with multiprocessing.Pool(processes=min(num_per_run, len(tasks))) as pool:
            pool.map(_process_sample, tasks)

        _collect_results(tags, output_path, outformat, file_detail, file_err, file_indel)

        for sample_id in samples:
            sample_file = f"{output_path}/{sample_id}"
            if os.path.exists(sample_file):
                os.remove(sample_file)

    else:
        # Multiple comma-separated input files
        tags = [Path(f).stem.split('.')[0] or Path(f).stem for f in input_files]
        tasks = [(f, output_path, tag) + common for f, tag in zip(input_files, tags)]

        with multiprocessing.Pool(processes=min(num_per_run, len(tasks))) as pool:
            pool.map(_process_sample, tasks)

        _collect_results(tags, output_path, outformat, file_detail, file_err, file_indel)


def main():
    parser = argparse.ArgumentParser(
        description='FunSeq2: Annotate and prioritize cancer somatic mutations'
    )
    parser.add_argument('infile', help='Input variants file (comma-separated for multiple)')
    parser.add_argument('maf', type=float, help='Minor allele frequency threshold')
    parser.add_argument('genome_mode', type=int, choices=[1, 2],
                        help='1: somatic; 2: germline | personal genome')
    parser.add_argument('informat', choices=['bed', 'vcf'], help='Input format [bed | vcf]')
    parser.add_argument('outformat', choices=['bed', 'vcf'], help='Output format [bed | vcf]')
    parser.add_argument('nc_mode', type=int, help='Only do non-coding if nc_mode=1')
    parser.add_argument('output_path', help='Output path')
    parser.add_argument('num_per_run', type=int, help='Number of parallel processes per run')
    parser.add_argument('cancer_type', help='Cancer type retrieved')
    parser.add_argument('score_cut', type=float, help='Non-coding candidate score cut')
    parser.add_argument('weight_mode', type=int, choices=[0, 1],
                        help='0: unweighted; 1: weighted scoring scheme')
    parser.add_argument('user_anno_dir', help='Directory for user-specific annotations')
    parser.add_argument('recur_db_use', help='Recurrence database use flag')
    parser.add_argument('sv_length_cut', help='SV length cut (integer or "inf")')
    parser.add_argument('gene_list', nargs='?', default='', help='Gene list file (optional)')
    parser.add_argument('--expression', default='', help='Expression file (optional)')
    parser.add_argument('--class', dest='class_file', default='', help='Sample class file (optional)')
    parser.add_argument('--exp_format', default='', help='Expression format: rpkm/raw (optional)')
    parser.add_argument('--config', default='config.txt', help='Configuration file path')

    args = parser.parse_args()

    exp_format = ''
    if args.expression and args.class_file:
        exp_format = args.exp_format if args.exp_format else 'rpkm'

    validate_parameters(args.maf, args.genome_mode, args.informat, args.outformat,
                        args.sv_length_cut, exp_format)

    os.makedirs(args.output_path, exist_ok=True)

    config_file = os.environ.get('FunSeqConfig', args.config)
    print(f"FunSeq is starting!\n############################\nConfig file:{config_file}\n")

    file_path, variables = read_config(config_file)

    # Detect GENCODE version from the gencode directory
    gencode_dir = variables.get('gencode', '')
    if gencode_dir:
        result = subprocess.run(
            f"ls {shlex.quote(gencode_dir)}/*.promoter.bed",
            shell=True, capture_output=True, text=True
        )
        if result.returncode == 0 and result.stdout.strip():
            gencode_file = result.stdout.strip().split('\n')[-1]
            parts = Path(gencode_file).stem.split('.')
            if len(parts) >= 2:
                gencode_v = parts[1]
                variables['cds'] = f"{gencode_dir}/gencode.{gencode_v}.cds.bed"
                variables['promoter'] = f"{gencode_dir}/gencode.{gencode_v}.promoter.bed"
                variables['intron'] = f"{gencode_dir}/gencode.{gencode_v}.intron.bed"
                variables['utr'] = f"{gencode_dir}/gencode.{gencode_v}.utr.bed"
                variables['coding_interval'] = f"{gencode_dir}/gencode.{gencode_v}.cds.interval"
                variables['coding_fasta'] = f"{gencode_dir}/gencode.{gencode_v}.cds.fa"

    # Validate required data files
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

    reference_file = variables.get('reference_file')
    if reference_file and not file_nonempty(reference_file):
        print(f"Warning: reference_file ({reference_file}) not found. Motif gain analysis will be skipped.\n")
        variables['reference_file'] = None

    if args.genome_mode == 2:
        check_file(variables.get('ancestral_file', ''), 'ancestral_file')

    if args.nc_mode == 0:
        check_file(variables.get('coding_interval', ''), 'coding_interval')
        check_file(variables.get('coding_fasta', ''), 'coding_fasta')

    file_detail = f"{args.output_path}/Output.{args.outformat}"
    file_recur = f"{args.output_path}/Recur.Summary"
    file_driver = f"{args.output_path}/Candidates.Summary"
    file_err = f"{args.output_path}/Error.log"
    file_indel = f"{args.output_path}/Output.indel.{args.outformat}"

    for f in [file_detail, file_recur, file_driver, file_err, file_indel]:
        if os.path.exists(f):
            os.remove(f)

    if args.expression and args.class_file:
        print("Differential Gene Expression Analysis ...\n")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script = os.path.join(os.path.dirname(script_dir), 'differential_gene_expression.r')
        subprocess.run(['Rscript', r_script, args.expression, args.class_file,
                        exp_format, args.output_path, file_err])

    input_files = args.infile.split(',')
    motif_p_value_cut = 4e-8

    process_files(
        input_files, args.output_path, variables,
        args.maf, args.informat, args.outformat, args.genome_mode,
        args.nc_mode, args.sv_length_cut, motif_p_value_cut,
        args.user_anno_dir, args.weight_mode, variables.get('weight_file', ''),
        args.gene_list, args.expression, args.class_file, exp_format,
        args.num_per_run, file_detail, file_err, file_indel
    )

    if file_nonempty(file_detail):
        Funseq_SNV.recur(file_detail, file_recur, file_driver, args.outformat,
                         variables.get('cancer_dir', ''), args.cancer_type, args.score_cut,
                         args.user_anno_dir, variables.get('weight_file', ''), args.recur_db_use)

    if file_nonempty(file_indel):
        Funseq_Indel.combine(file_indel, args.outformat)


if __name__ == '__main__':
    main()
