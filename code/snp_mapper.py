#!/usr/bin/env python3
"""Python replacement for VAT's snpMapper tool.

Usage: snp_mapper.py <annotation.interval> <annotation.fa>
Reads VCF from stdin, writes annotated VCF to stdout.
"""

import sys
from collections import defaultdict

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def load_intervals(interval_file):
    """
    Returns:
      chrom_index: chrom -> list of transcript dicts
      gene_tx_count: gene_id -> total number of transcripts
    """
    transcripts = []
    with open(interval_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            name = fields[0]
            chrom = fields[1]
            strand = fields[2]
            exon_starts = [int(x) for x in fields[6].rstrip(',').split(',') if x]
            exon_ends   = [int(x) for x in fields[7].rstrip(',').split(',') if x]
            parts = name.split('|')
            transcripts.append({
                'name':            name,
                'gene_id':         parts[0] if len(parts) > 0 else name,
                'transcript_id':   parts[1] if len(parts) > 1 else '',
                'gene_name':       parts[2] if len(parts) > 2 else name,
                'transcript_name': parts[3] if len(parts) > 3 else '',
                'chrom':           chrom,
                'strand':          strand,
                'exon_starts':     exon_starts,
                'exon_ends':       exon_ends,
            })

    chrom_index = defaultdict(list)
    gene_tx_count = defaultdict(int)
    for t in transcripts:
        chrom_index[t['chrom']].append(t)
        gene_tx_count[t['gene_id']] += 1

    return chrom_index, gene_tx_count


def load_fasta(fasta_file):
    """Returns dict: full header string -> CDS sequence."""
    seqs = {}
    name = None
    buf = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def fasta_key(t):
    """Build the FASTA header key that matches load_fasta keys."""
    parts = [t['gene_id'], t['transcript_id'], t['gene_name'], t['transcript_name'],
             t['chrom'], t['strand']]
    for es, ee in zip(t['exon_starts'], t['exon_ends']):
        parts.append(str(es))
        parts.append(str(ee))
    return '|'.join(parts)


def get_tx_coord(t, genomic_pos):
    """
    Map 0-based genomic_pos to 0-based CDS transcript coordinate.
    Returns None if pos is not in any exon.
    """
    offset = 0
    for es, ee in zip(t['exon_starts'], t['exon_ends']):
        if es <= genomic_pos < ee:
            offset += (genomic_pos - es)
            if t['strand'] == '-':
                total = sum(e - s for s, e in zip(t['exon_starts'], t['exon_ends']))
                return total - offset - 1
            return offset
        offset += (ee - es)
    return None


TYPE_PRIORITY = {'prematureStop': 3, 'nonsynonymous': 2, 'removedStop': 1, 'synonymous': 0}


def annotate_variant(chrom, pos_1based, ref, alt, chrom_index, gene_tx_count, seqs):
    """
    Returns VA annotation string for a single SNP, or None.
    pos_1based: VCF 1-based position.
    """
    if len(ref) != 1 or len(alt) != 1:
        return None

    pos = pos_1based - 1  # 0-based

    # Find all transcripts where pos falls inside an exon
    overlapping = [
        t for t in chrom_index.get(chrom, [])
        if any(es <= pos < ee for es, ee in zip(t['exon_starts'], t['exon_ends']))
    ]
    if not overlapping:
        return None

    # Group by gene
    by_gene = defaultdict(list)
    for t in overlapping:
        by_gene[t['gene_id']].append(t)

    va_parts = []
    for gene_id in sorted(by_gene):
        txs = by_gene[gene_id]
        gene_name = txs[0]['gene_name']
        strand     = txs[0]['strand']

        # For minus-strand, the CDS sequence is already 5'->3' (revcomped),
        # so we need the complement of the genomic alleles.
        if strand == '-':
            cds_ref = revcomp(ref)
            cds_alt = revcomp(alt)
        else:
            cds_ref = ref
            cds_alt = alt

        tx_details = []
        for t in sorted(txs, key=lambda x: x['transcript_name']):
            tx_coord = get_tx_coord(t, pos)
            if tx_coord is None:
                continue
            raw_seq = seqs.get(fasta_key(t))
            if not raw_seq or tx_coord >= len(raw_seq):
                continue
            # The FASTA stores the plus-strand sequence concatenated in genomic
            # (ascending) exon order.  For minus-strand genes the actual CDS
            # (5'->3') is therefore the reverse complement of the stored sequence.
            cds = revcomp(raw_seq) if t['strand'] == '-' else raw_seq
            if cds[tx_coord].upper() != cds_ref.upper():
                continue  # reference mismatch — skip this transcript

            alt_seq = cds[:tx_coord] + cds_alt + cds[tx_coord+1:]
            prot_ref = translate(cds)
            prot_alt = translate(alt_seq)

            stops_ref = prot_ref.count('*')
            stops_alt = prot_alt.count('*')
            if   stops_alt > stops_ref:  mut = 'prematureStop'
            elif stops_alt < stops_ref:  mut = 'removedStop'
            elif prot_ref == prot_alt:   mut = 'synonymous'
            else:                        mut = 'nonsynonymous'

            cds_len = len(cds)
            nt_pos  = tx_coord + 1               # 1-based
            aa_pos  = tx_coord // 3 + 1          # 1-based
            ref_aa  = prot_ref[aa_pos-1] if aa_pos-1 < len(prot_ref) else 'X'
            alt_aa  = prot_alt[aa_pos-1] if aa_pos-1 < len(prot_alt) else 'X'

            tx_details.append((mut,
                                t['transcript_name'], t['transcript_id'],
                                cds_len, nt_pos, aa_pos, ref_aa, alt_aa))

        if not tx_details:
            continue

        dominant = max(tx_details, key=lambda x: TYPE_PRIORITY.get(x[0], 0))[0]
        n_affected = len(tx_details)
        total_tx   = gene_tx_count.get(gene_id, n_affected)

        # Build transcript portion: name:id:len_ntpos_aapos_refAA->altAA:...
        tx_str = ':'.join(
            f"{tx_name}:{tx_id}:{cds_len}_{nt_pos}_{aa_pos}_{ref_aa}->{alt_aa}"
            for _, tx_name, tx_id, cds_len, nt_pos, aa_pos, ref_aa, alt_aa in tx_details
        )

        va_parts.append(
            f"1:{gene_name}:{gene_id}:{strand}:{dominant}:{n_affected}/{total_tx}:{tx_str}"
        )

    return ','.join(va_parts) if va_parts else None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <annotation.interval> <annotation.fa>", file=sys.stderr)
        sys.exit(1)

    print("Loading intervals...", file=sys.stderr)
    chrom_index, gene_tx_count = load_intervals(sys.argv[1])
    print("Loading CDS sequences...", file=sys.stderr)
    seqs = load_fasta(sys.argv[2])
    print(f"Loaded {len(seqs)} CDS sequences", file=sys.stderr)

    header_printed = False
    pending_comments = []

    for raw in sys.stdin:
        line = raw.rstrip('\n')
        if line.startswith('##'):
            pending_comments.append(line)
            continue
        if line.startswith('#CHROM'):
            for c in pending_comments:
                print(c)
            print('##INFO=<ID=VA,Number=.,Type=String,Description="Variant Annotation">')
            print(line)
            header_printed = True
            continue

        if not header_printed:
            for c in pending_comments:
                print(c)
            print('##INFO=<ID=VA,Number=.,Type=String,Description="Variant Annotation">')
            header_printed = True

        fields = line.split('\t')
        if len(fields) < 5:
            continue
        chrom, pos_s, vcf_id, ref, alt_field = fields[0], fields[1], fields[2], fields[3], fields[4]
        qual   = fields[5] if len(fields) > 5 else '.'
        filt   = fields[6] if len(fields) > 6 else '.'
        info   = fields[7] if len(fields) > 7 else '.'

        try:
            pos = int(pos_s)
        except ValueError:
            continue

        va_annotations = []
        for alt in alt_field.split(','):
            va = annotate_variant(chrom, pos, ref, alt, chrom_index, gene_tx_count, seqs)
            if va:
                va_annotations.append(va)

        if va_annotations:
            va_str  = ','.join(va_annotations)
            new_info = f"{info};VA={va_str}" if info not in ('.', '') else f"VA={va_str}"
            out = [chrom, pos_s, vcf_id, ref, alt_field, qual, filt, new_info]
            if len(fields) > 8:
                out.extend(fields[8:])
            print('\t'.join(out))


if __name__ == '__main__':
    main()
