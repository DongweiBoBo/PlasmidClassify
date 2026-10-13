#!/usr/bin/env python3
"""
Extract plasmid contigs from Flye assembly results, separated by circularity.

Reads 3CAC classification to identify plasmids, assembly_info.txt to determine
circularity, and extracts sequences from the assembly FASTA file.

Usage:
    python extract_plasmids_flye.py \
        --classification Sample.flye/3CAC_classification.fasta \
        --assembly_info assembly_info.txt \
        --assembly_fasta assembly.fasta \
        --output_dir Sample.flye
"""

import argparse
import os
import sys


def read_classification(filepath):
    """Read 3CAC classification file. Returns dict: contig_name -> category"""
    classification = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.rsplit(',', 1)
            if len(parts) == 2:
                name, category = parts
                classification[name.strip()] = category.strip()
    return classification


def read_circular_info(filepath):
    """Read Flye assembly_info.txt to get circularity info.
    Returns dict: contig_name -> is_circular (bool)
    The circ. column (4th column, 0-indexed: 3) contains 'Y' or 'N'.
    """
    circular = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                name = parts[0]
                is_circ = parts[3].strip().upper() == 'Y'
                circular[name] = is_circ
    return circular


def read_fasta(filepath):
    """Read FASTA file. Returns dict: contig_name -> sequence"""
    sequences = {}
    name = None
    seq_parts = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(seq_parts)
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            sequences[name] = ''.join(seq_parts)

    return sequences


def write_fasta(filepath, names, sequences):
    """Write selected sequences to a FASTA file."""
    with open(filepath, 'w') as f:
        for name in names:
            if name in sequences:
                seq = sequences[name]
                f.write(f'>{name}\n')
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Extract plasmid contigs separated by circularity')
    parser.add_argument('--classification', required=True,
                        help='3CAC classification result file')
    parser.add_argument('--assembly_info', required=True,
                        help='Flye assembly_info.txt')
    parser.add_argument('--assembly_fasta', required=True,
                        help='Original assembly FASTA file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Read classification results
    classification = read_classification(args.classification)
    plasmid_contigs = {name for name, cat in classification.items() if cat == 'plasmid'}
    print(f"Total contigs: {len(classification)}")
    print(f"Plasmid contigs: {len(plasmid_contigs)}")

    # 2. Read circularity info
    circular_info = read_circular_info(args.assembly_info)

    # 3. Separate plasmids by circularity
    circular_plasmids = []
    linear_plasmids = []

    for contig in sorted(plasmid_contigs, key=lambda x: int(x.split('_')[1])):
        is_circ = circular_info.get(contig, False)
        if is_circ:
            circular_plasmids.append(contig)
        else:
            linear_plasmids.append(contig)

    print(f"Circular plasmids: {len(circular_plasmids)}")
    print(f"Non-circular plasmids: {len(linear_plasmids)}")

    # 4. Read assembly sequences
    print(f"Reading assembly FASTA...")
    sequences = read_fasta(args.assembly_fasta)
    print(f"Loaded {len(sequences)} sequences")

    # 5. Write output files
    circular_file = os.path.join(args.output_dir, 'plasmid_circular.fasta')
    linear_file = os.path.join(args.output_dir, 'plasmid_linear.fasta')

    write_fasta(circular_file, circular_plasmids, sequences)
    write_fasta(linear_file, linear_plasmids, sequences)

    print(f"\nCircular plasmids written to: {circular_file}")
    print(f"Non-circular plasmids written to: {linear_file}")

    # 6. Summary statistics
    if circular_plasmids:
        circ_lengths = [len(sequences.get(c, '')) for c in circular_plasmids]
        print(f"\nCircular plasmids - count: {len(circ_lengths)}, "
              f"length range: {min(circ_lengths)}-{max(circ_lengths)} bp, "
              f"total: {sum(circ_lengths)} bp")
    if linear_plasmids:
        lin_lengths = [len(sequences.get(c, '')) for c in linear_plasmids]
        print(f"Non-circular plasmids - count: {len(lin_lengths)}, "
              f"length range: {min(lin_lengths)}-{max(lin_lengths)} bp, "
              f"total: {sum(lin_lengths)} bp")


if __name__ == '__main__':
    main()
