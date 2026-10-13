#!/usr/bin/env python3
"""
Generate scaffolds.paths file for 3CAC from SPAdes output.
"""
import sys

gfa_file = sys.argv[1]       
fasta_file = sys.argv[2]    
output_file = sys.argv[3]   

edge_ids = []
fasta_headers = []

with open(gfa_file) as gf:
    for line in gf:
        if line.startswith('S'):
            edge_ids.append(line.strip().split('\t')[1])
        else:
            break 

with open(fasta_file) as ff:
    for line in ff:
        if line.startswith('>'):
            fasta_headers.append(line.strip()[1:]) 

with open(output_file, 'w') as out:
    for i in range(min(len(edge_ids), len(fasta_headers))):
        header = fasta_headers[i]
        edge_id = edge_ids[i]
        out.write(f"{header}\n")
        out.write(f"{edge_id}+\n")
        out.write(f"{header}'\n")
        out.write(f"{edge_id}-\n")
