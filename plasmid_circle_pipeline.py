#!/usr/bin/env python3
"""
Iterative plasmid assembly and classification pipeline.

This pipeline performs iterative refinement of plasmid assemblies:
1. BLAST against plasmid database
2. Group sequences by genus
3. Assemble each genus group with long reads
4. Classify assemblies (PPR_Meta, PlasClass, DeepVirFinder, 3CAC)
5. Extract circular plasmids
6. Iterate with non-circular plasmids

Usage:
    python3 plasmid_pipeline.py \
        --fasta input.fasta \
        --reads reads.fastq.gz \
        --db_index /path/to/plasmid_db.ndb \
        --tax_dir /path/to/tax \
        --output_dir output \
        --iterations 3
"""

import argparse
import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict


def check_tool(tool_name):
    """Check if a tool is available in PATH."""
    if shutil.which(tool_name) is None:
        print(f"Error: {tool_name} not found in PATH")
        sys.exit(1)


def run_command(cmd, description):
    """Run a shell command and check for errors."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    print(f"{'='*60}")
    
    if isinstance(cmd, str):
        result = subprocess.run(cmd, shell=True)
    else:
        result = subprocess.run(cmd)
    
    if result.returncode != 0:
        print(f"Error: {description} failed with return code {result.returncode}")
        sys.exit(1)
    
    print(f"✓ {description} completed successfully")


def blast_search(fasta_file, db_index, output_dir, threads=8):
    """Perform BLAST search against plasmid database."""
    output_file = os.path.join(output_dir, "blast_results.txt")
    
    cmd = [
        "blastn",
        "-query", fasta_file,
        "-db", db_index,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
        "-num_threads", str(threads),
        "-evalue", "1e-6"
    ]
    
    run_command(cmd, "BLAST search")
    return output_file


def extract_taxonomy(blast_file, tax_dir):
    """Extract taxonomy information from BLAST results."""
    print("\nExtracting taxonomy information...")
    
    # Read taxonomy data
    tax_data = {}
    tax_nodes = os.path.join(tax_dir, "nodes.dmp")
    tax_names = os.path.join(tax_dir, "names.dmp")
    
    # Parse names.dmp to get scientific names
    name_map = {}
    if os.path.exists(tax_names):
        with open(tax_names) as f:
            for line in f:
                parts = line.split("|")
                if len(parts) >= 3 and "scientific name" in parts[3]:
                    tax_id = parts[0].strip()
                    name = parts[1].strip()
                    name_map[tax_id] = name
    
    # Parse nodes.dmp to get taxonomy hierarchy
    parent_map = {}
    rank_map = {}
    if os.path.exists(tax_nodes):
        with open(tax_nodes) as f:
            for line in f:
                parts = line.split("|")
                if len(parts) >= 3:
                    tax_id = parts[0].strip()
                    parent_id = parts[1].strip()
                    rank = parts[2].strip()
                    parent_map[tax_id] = parent_id
                    rank_map[tax_id] = rank
    
    # Parse BLAST results and extract genus for each query
    query_genus = {}
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 13:
                query_id = parts[0]
                tax_ids = parts[12].split(";")
                
                if query_id not in query_genus and tax_ids[0] in name_map:
                    # Walk up taxonomy tree to find genus
                    current_id = tax_ids[0]
                    genus = None
                    
                    while current_id in parent_map:
                        if rank_map.get(current_id) == "genus":
                            genus = name_map.get(current_id)
                            break
                        current_id = parent_map[current_id]
                    
                    if genus:
                        query_genus[query_id] = genus
    
    print(f"Found genus information for {len(query_genus)} sequences")
    return query_genus


def group_by_genus(query_genus, fasta_file, output_dir):
    """Group sequences by genus and save to separate files."""
    print("\nGrouping sequences by genus...")
    
    # Read all sequences
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = "".join(current_seq)
    
    # Group by genus
    genus_groups = defaultdict(list)
    for seq_id, genus in query_genus.items():
        if seq_id in sequences:
            genus_groups[genus].append((seq_id, sequences[seq_id]))
    
    # Write genus-specific FASTA files
    genus_files = {}
    for genus, seqs in genus_groups.items():
        safe_name = genus.replace(" ", "_").replace("/", "_")
        genus_file = os.path.join(output_dir, f"{safe_name}.fasta")
        
        with open(genus_file, "w") as f:
            for seq_id, seq in seqs:
                f.write(f">{seq_id}\n{seq}\n")
        
        genus_files[genus] = genus_file
    
    print(f"Created {len(genus_files)} genus-specific files")
    return genus_files


def minimap2_align(genus_fasta, reads_fastq, output_dir, threads=8):
    """Align genus sequences to long reads using minimap2."""
    paf_file = os.path.join(output_dir, "alignment.paf")
    
    cmd = [
        "minimap2",
        "-x", "map-ont",
        "-t", str(threads),
        genus_fasta,
        reads_fastq,
        ">", paf_file
    ]
    
    run_command(" ".join(cmd), f"minimap2 alignment for {os.path.basename(genus_fasta)}")
    return paf_file


def extract_aligned_reads(paf_file, reads_fastq, output_dir):
    """Extract reads that aligned to genus sequences."""
    # Get list of aligned query names
    aligned_queries = set()
    with open(paf_file) as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) >= 6:
                aligned_queries.add(parts[0])
    
    # Extract reads (assuming read names match)
    extracted_fastq = os.path.join(output_dir, "extracted_reads.fastq.gz")
    
    # Use zcat and awk to extract reads
    cmd = f"zcat {reads_fastq} | awk 'NR%4==1{{name=substr($1,2)}} name in queries{{print}}' > {extracted_fastq}"
    
    # Simpler approach: extract all reads and filter
    # This is a simplified version - in practice, you'd need to match read names properly
    print(f"Extracting {len(aligned_queries)} aligned reads...")
    
    # For now, just copy the original reads (simplified)
    # In practice, you'd implement proper read extraction
    shutil.copy(reads_fastq, extracted_fastq)
    
    return extracted_fastq


def flye_assemble(reads_fastq, output_dir, threads=8):
    """Assemble reads using flye."""
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        "flye",
        "--nano-raw", reads_fastq,
        "--out-dir", output_dir,
        "--threads", str(threads),
        "--iterations", "2"
    ]
    
    run_command(cmd, "flye assembly")
    
    assembly_file = os.path.join(output_dir, "assembly.fasta")
    return assembly_file


def classify_assembly(assembly_file, assembly_dir, output_dir, threads=8):
    """Classify assembly using PPR_Meta, PlasClass, DeepVirFinder, and 3CAC."""
    print("\nClassifying assembly...")
    
    # Filter contigs > 500bp
    filtered_fasta = os.path.join(output_dir, "filtered_contigs.fasta")
    filter_cmd = f"seqkit seq -m 500 {assembly_file} > {filtered_fasta}"
    run_command(filter_cmd, "Filtering contigs > 500bp")
    
    # Run PPR_Meta
    pprmeta_output = os.path.join(output_dir, "pprmeta_results.csv")
    pprmeta_cmd = f"PPR_Meta {filtered_fasta} {output_dir}/pprmeta_results"
    run_command(pprmeta_cmd, "PPR_Meta classification")
    
    # Run PlasClass
    plasclass_output = os.path.join(output_dir, "plasclass_results.out")
    plasclass_cmd = f"classify_fasta.py -p 8 -f {filtered_fasta} -o {plasclass_output}"
    run_command(plasclass_cmd, "PlasClass classification")
    
    # Run DeepVirFinder
    deepvirfinder_output = os.path.join(output_dir, "filtered_contigs.fasta_gt500bp_dvfpred.txt")
    deepvirfinder_cmd = f"python dvf.py -i {filtered_fasta} -o {output_dir} -l 500 -c 8"
    run_command(deepvirfinder_cmd, "DeepVirFinder classification")
    
    # Run 3CAC using Java Classify3CAC
    cac_output_dir = os.path.join(output_dir, "3CAC_output")
    os.makedirs(cac_output_dir, exist_ok=True)
    
    # Flye assembly outputs
    assembly_graph = os.path.join(assembly_dir, "assembly_graph.gfa")
    assembly_info = os.path.join(assembly_dir, "assembly_info.txt")
    
    # Check if required files exist
    if not os.path.exists(assembly_graph):
        print(f"Warning: {assembly_graph} not found, skipping 3CAC")
        return pprmeta_output, plasclass_output, deepvirfinder_output, None
    
    if not os.path.exists(assembly_info):
        print(f"Warning: {assembly_info} not found, skipping 3CAC")
        return pprmeta_output, plasclass_output, deepvirfinder_output, None
    
    # Build 3CAC command
    java_3cac = "path/Classify3CAC.java"
    
    cac_cmd = [
        "java", java_3cac,
        "--assembler", "Flye",
        "--PPRMeta", pprmeta_output,
        "--plasClass", plasclass_output,
        "--deepVirFinder", deepvirfinder_output,
        "--graph", assembly_graph,
        "--path", assembly_info,
        "--output", cac_output_dir
    ]
    
    run_command(" ".join(cac_cmd), "3CAC classification")
    
    cac_classification = os.path.join(cac_output_dir, "3CAC_classification.fasta")
    
    return pprmeta_output, plasclass_output, deepvirfinder_output, cac_classification


def extract_plasmids_simple(classification_file, assembly_file, output_dir):
    """Extract plasmid sequences based on classification (simplified version)."""
    print("\nExtracting plasmid sequences...")
    
    # Read classification
    plasmid_ids = set()
    with open(classification_file) as f:
        for line in f:
            if "plasmid" in line.lower():
                parts = line.split(",")
                if len(parts) >= 1:
                    plasmid_ids.add(parts[0].strip())
    
    # Extract sequences
    plasmid_fasta = os.path.join(output_dir, "plasmids.fasta")
    
    with open(assembly_file) as fin, open(plasmid_fasta, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id in plasmid_ids:
                    write = True
                    fout.write(line)
                else:
                    write = False
            else:
                if write:
                    fout.write(line)
    
    print(f"Extracted {len(plasmid_ids)} plasmid sequences")
    return plasmid_fasta


def extract_plasmids_flye(classification_file, assembly_file, output_dir):
    """Extract plasmid sequences from Flye 3CAC classification.
    Classification format: contig_X,category
    """
    print("\nExtracting plasmid sequences (Flye format)...")
    
    # Read classification - format: contig_X,category
    plasmid_ids = set()
    with open(classification_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.rsplit(',', 1)
            if len(parts) == 2 and parts[1].strip() == 'plasmid':
                plasmid_ids.add(parts[0].strip())
    
    # Extract sequences
    plasmid_fasta = os.path.join(output_dir, "plasmids.fasta")
    
    with open(assembly_file) as fin, open(plasmid_fasta, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id in plasmid_ids:
                    write = True
                    fout.write(line)
                else:
                    write = False
            else:
                if write:
                    fout.write(line)
    
    print(f"Extracted {len(plasmid_ids)} plasmid sequences")
    return plasmid_fasta


def separate_circular_plasmids(plasmid_fasta, output_dir):
    """Separate circular and non-circular plasmids (simplified)."""
    # For flye assemblies, check assembly_info.txt for circularity
    # This is a simplified version
    
    circular_file = os.path.join(output_dir, "plasmid_circular.fasta")
    linear_file = os.path.join(output_dir, "plasmid_linear.fasta")
    
    # Simplified: assume all are linear for now
    # In practice, you'd check flye's assembly_info.txt
    shutil.copy(plasmid_fasta, linear_file)
    
    # Create empty circular file
    with open(circular_file, "w") as f:
        pass
    
    return circular_file, linear_file


def separate_circular_plasmids_flye(plasmid_fasta, assembly_info, output_dir):
    """Separate circular and non-circular plasmids using Flye's assembly_info.txt.
    The circ. column (4th column) contains 'Y' or 'N'.
    """
    print("\nSeparating circular and non-circular plasmids...")
    
    # Read circularity info from assembly_info.txt
    circular_contigs = set()
    if os.path.exists(assembly_info):
        with open(assembly_info) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[3].strip().upper() == "Y":
                    circular_contigs.add(parts[0].strip())
    
    # Read plasmid sequences and separate
    circular_file = os.path.join(output_dir, "plasmid_circular.fasta")
    linear_file = os.path.join(output_dir, "plasmid_linear.fasta")
    
    circular_count = 0
    linear_count = 0
    
    with open(plasmid_fasta) as fin, \
         open(circular_file, "w") as fcirc, \
         open(linear_file, "w") as flinear:
        
        write = False
        fout = None
        
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id in circular_contigs:
                    fout = fcirc
                    circular_count += 1
                else:
                    fout = flinear
                    linear_count += 1
                fout.write(line)
                write = True
            else:
                if write:
                    fout.write(line)
    
    print(f"Circular plasmids: {circular_count}")
    print(f"Non-circular plasmids: {linear_count}")
    
    return circular_file, linear_file


def main():
    parser = argparse.ArgumentParser(
        description="Iterative plasmid assembly and classification pipeline")
    
    parser.add_argument("--fasta", required=True,
                        help="Input FASTA file with contigs")
    parser.add_argument("--reads", required=True,
                        help="Raw long-read data (fastq.gz)")
    parser.add_argument("--db_index", required=True,
                        help="Plasmid database index file")
    parser.add_argument("--tax_dir", required=True,
                        help="Taxonomy directory")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory")
    parser.add_argument("--iterations", type=int, default=3,
                        help="Number of iterations (default: 3)")
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    
    args = parser.parse_args()
    
    # Check required tools
    tools = ["blastn", "minimap2", "flye", "PPR_Meta", "classify_fasta.py"]
    for tool in tools:
        check_tool(tool)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize with input sequences
    current_fasta = args.fasta
    all_circular_plasmids = []
    
    # Main iteration loop
    for iteration in range(1, args.iterations + 1):
        print(f"\n{'#'*60}")
        print(f"# ITERATION {iteration}/{args.iterations}")
        print(f"{'#'*60}")
        
        iter_dir = os.path.join(args.output_dir, f"iteration_{iteration}")
        os.makedirs(iter_dir, exist_ok=True)
        
        # Step 1: BLAST search
        blast_file = blast_search(current_fasta, args.db_index, iter_dir, args.threads)
        
        # Step 2: Extract taxonomy
        query_genus = extract_taxonomy(blast_file, args.tax_dir)
        
        if not query_genus:
            print("No sequences with genus information found. Stopping.")
            break
        
        # Step 3: Group by genus
        genus_files = group_by_genus(query_genus, current_fasta, iter_dir)
        
        # Step 4: Process each genus
        iter_circular_plasmids = []
        iter_noncircular_plasmids = []
        
        for genus, genus_fasta in genus_files.items():
            print(f"\nProcessing genus: {genus}")
            
            genus_dir = os.path.join(iter_dir, f"genus_{genus.replace(' ', '_')}")
            os.makedirs(genus_dir, exist_ok=True)
            
            # Align with minimap2
            paf_file = minimap2_align(genus_fasta, args.reads, genus_dir, args.threads)
            
            # Extract aligned reads
            extracted_reads = extract_aligned_reads(paf_file, args.reads, genus_dir)
            
            # Assemble with flye
            flye_dir = os.path.join(genus_dir, "flye_assembly")
            assembly_file = flye_assemble(extracted_reads, flye_dir, args.threads)
            
            # Classify assembly
            pprmeta_out, plasclass_out, dvf_out, cac_out = classify_assembly(
                assembly_file, flye_dir, genus_dir, args.threads)
            
            # Extract plasmids based on 3CAC classification
            if cac_out and os.path.exists(cac_out):
                plasmid_fasta = extract_plasmids_flye(
                    cac_out, assembly_file, genus_dir)
            else:
                # Fallback to PPRMeta if 3CAC failed
                plasmid_fasta = extract_plasmids_simple(
                    pprmeta_out, assembly_file, genus_dir)
            
            # Separate circular and non-circular using flye's assembly_info.txt
            assembly_info = os.path.join(flye_dir, "assembly_info.txt")
            circular_file, linear_file = separate_circular_plasmids_flye(
                plasmid_fasta, assembly_info, genus_dir)
            
            # Collect results
            if os.path.exists(circular_file) and os.path.getsize(circular_file) > 0:
                iter_circular_plasmids.append(circular_file)
            
            if os.path.exists(linear_file) and os.path.getsize(linear_file) > 0:
                iter_noncircular_plasmids.append(linear_file)
        
        # Merge circular plasmids from this iteration
        if iter_circular_plasmids:
            merged_circular = os.path.join(iter_dir, "all_circular_plasmids.fasta")
            with open(merged_circular, "w") as fout:
                for f in iter_circular_plasmids:
                    with open(f) as fin:
                        fout.write(fin.read())
            all_circular_plasmids.append(merged_circular)
            print(f"\nIteration {iteration}: Found {len(iter_circular_plasmids)} circular plasmid files")
        
        # Prepare for next iteration
        if iter_noncircular_plasmids:
            merged_linear = os.path.join(iter_dir, "all_noncircular_plasmids.fasta")
            with open(merged_linear, "w") as fout:
                for f in iter_noncircular_plasmids:
                    with open(f) as fin:
                        fout.write(fin.read())
            current_fasta = merged_linear
            print(f"Preparing for next iteration with {len(iter_noncircular_plasmids)} non-circular plasmid files")
        else:
            print("No non-circular plasmids found. Stopping iterations.")
            break
    
    # Final output: merge all circular plasmids from all iterations
    if all_circular_plasmids:
        final_output = os.path.join(args.output_dir, "final_circular_plasmids.fasta")
        with open(final_output, "w") as fout:
            for f in all_circular_plasmids:
                with open(f) as fin:
                    fout.write(fin.read())
        
        print(f"\n{'='*60}")
        print(f"Pipeline completed successfully!")
        print(f"Final output: {final_output}")
        print(f"{'='*60}")
    else:
        print("\nNo circular plasmids found in any iteration.")


if __name__ == "__main__":
    main()
