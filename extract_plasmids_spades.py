#!/usr/bin/env python3
"""
Extract plasmid contigs from SPAdes 3CAC classification results.

Usage:
    python3 extract_plasmids_spades.py \
        --classification Sample.spades/3CAC_classification.fasta \
        --assembly_fasta Sample.assembly.fasta \
        --output Sample.spades/plasmid_linear.fasta
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description='Extract SPAdes plasmid contigs')
    parser.add_argument('--classification', required=True, help='3CAC classification file')
    parser.add_argument('--assembly_fasta', required=True, help='Assembly FASTA file')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    args = parser.parse_args()

    # 1. Read classification, collect plasmid contig IDs (e.g., NODE_1, NODE_2, ...)
    plasmid_ids = set()
    with open(args.classification) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, category = line.rsplit(',', 1)
            if category.strip() == 'plasmid':
                plasmid_ids.add(name.strip())

    print(f"Plasmid contigs in classification: {len(plasmid_ids)}")

    # 2. Read FASTA and extract plasmid sequences
    #    FASTA headers are like >NODE_1_length_1135_cov_0.459325
    #    Match by extracting NODE_X prefix
    count = 0
    with open(args.assembly_fasta) as fin, open(args.output, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                header = line.strip()
                # Extract NODE_X from >NODE_X_length_Y_cov_Z
                contig_id = header[1:].split('_')[0] + '_' + header[1:].split('_')[1]
                # More robust: take everything before _length_
                parts = header[1:].split('_length_')
                contig_id = parts[0]  # NODE_X

                if contig_id in plasmid_ids:
                    write = True
                    fout.write(header + '\n')
                    count += 1
                else:
                    write = False
            else:
                if write:
                    fout.write(line)

    print(f"Extracted {count} plasmid sequences to: {args.output}")


if __name__ == '__main__':
    main()
