#!/usr/bin/env python3
"""
3CAC-like classifier for Unicycler assembly results.

Usage:
    python3 classify_3CAC_unicycler.py \
        --PPRMeta <PPRMeta_csv> \
        --plasClass <PlasClass_output> \
        --deepVirFinder <DeepVirFinder_output> \
        --graph <assembly.gfa> \
        --output <output_directory>
"""

import sys
import os
import argparse
from collections import defaultdict


def extract_contig_id(name):
    """Extract numeric contig ID from various name formats.
    E.g., '1 length=564934 depth=1.45x' -> 1, '29' -> 29
    """
    return int(name.strip().split()[0])


def read_plasclass(filepath):
    """Read PlasClass predictions. Returns dict: contig_id -> score
    PlasClass output may be tab or space separated: contig_id\tscore
    """
    plas = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                parts = line.split()
            cid = extract_contig_id(parts[0])
            score = float(parts[1])
            plas[cid] = score
    return plas


def read_deepvirfinder(filepath):
    """Read DeepVirFinder predictions. Returns dict: contig_id -> (score, pvalue)
    DVF output is tab-separated: name\tlen\tscore\tpvalue
    The name field may contain spaces (e.g., '29 length=132138 depth=0.33x')
    """
    dvf = {}
    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            cid = extract_contig_id(parts[0])
            score = float(parts[2])
            pvalue = float(parts[3])
            dvf[cid] = (score, pvalue)
    return dvf


def get_cate_pprmeta(fields):
    """Determine category from PPR-Meta classification.
    Returns: 1=phage, 2=plasmid, 3=chromosome, 4=uncertain
    """
    PPRMeta_score = 0.7
    source = fields[5].strip().lower()
    phage_score = float(fields[2])
    chrom_score = float(fields[3])
    plasmid_score = float(fields[4])

    if source == 'phage' and phage_score >= PPRMeta_score:
        return 1
    elif source == 'plasmid' and plasmid_score >= PPRMeta_score:
        return 2
    elif source == 'chromosome' and chrom_score >= PPRMeta_score:
        return 3
    else:
        return 4


def read_pprmeta(filepath):
    """Read PPR-Meta classification.
    Returns list of (contig_name, contig_id, pprmeta_category)
    """
    entries = []
    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split(',')
            name = fields[0].strip()
            cid = extract_contig_id(name)
            cate = get_cate_pprmeta(fields)
            entries.append((name, cid, cate))
    return entries


def initial_classification(entries, plas, dvf):
    """Perform initial classification using PPRMeta + PlasClass + DeepVirFinder.
    Same logic as 3CAC Java Classify3CAC.initialClassification()
    """
    initial_class = []
    num_phageC = 0
    num_plasC = 0
    num_uncertain = 0

    for name, cid, cate in entries:
        if cate == 1:  # PPRMeta says phage
            if cid in dvf:
                score, pvalue = dvf[cid]
                if pvalue <= 0.03:
                    cate = 1  # keep as phage
                elif score <= 0.5:
                    cate = 3  # reclassify as chromosome
                    num_phageC += 1
                else:
                    cate = 4  # uncertain
                    num_uncertain += 1
            else:
                cate = 4
                num_uncertain += 1

        elif cate == 2:  # PPRMeta says plasmid
            if cid in plas:
                score = plas[cid]
                if score >= 0.7:
                    cate = 2  # keep as plasmid
                elif score <= 0.3:
                    cate = 3  # reclassify as chromosome
                    num_plasC += 1
                else:
                    cate = 4  # uncertain
                    num_uncertain += 1
            else:
                cate = 4
                num_uncertain += 1

        # cate 3 (chromosome) and 4 (uncertain) remain unchanged
        initial_class.append((name, cid, cate))

    return initial_class, num_phageC, num_plasC, num_uncertain


def read_assembly_graph(gfa_file):
    """Read assembly graph from GFA file.
    Returns:
        edges_length: dict edge_id -> sequence_length
        edges_adj: dict edge_id -> set of adjacent edge_ids
        max_edge_id: maximum edge ID
    """
    edges_length = {}
    edges_adj = defaultdict(set)
    max_edge_id = 0

    with open(gfa_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                eid = int(parts[1])
                if eid > max_edge_id:
                    max_edge_id = eid
                edges_length[eid] = len(parts[2])
            elif parts[0] == 'L':
                eid1 = int(parts[1])
                eid2 = int(parts[3])
                if eid1 != eid2:
                    edges_adj[eid1].add(eid2)
                    edges_adj[eid2].add(eid1)
            else:
                if parts[0] not in ('S', 'L', ''):
                    break

    return edges_length, edges_adj, max_edge_id


def contig_category_to_edges(initial_class, edges_cate):
    """Assign contig categories to their corresponding edges.
    For Unicycler: contig N maps to edge N (one-to-one).
    """
    for name, cid, cate in initial_class:
        if cate != 4:
            if edges_cate[cid] == 4:
                edges_cate[cid] = cate
            elif edges_cate[cid] != cate:
                edges_cate[cid] = -1  # conflict

    # Set conflict edges to uncertain
    for eid in edges_cate:
        if edges_cate[eid] == -1:
            edges_cate[eid] = 4


def correction(edges_adj, edges_cate, max_edge_id):
    """Correction step: correct classified edges based on neighbor consensus.
    Same logic as 3CAC Java Classify3CAC.correction()
    """
    correct_edges = 0
    loop_num = 0

    correct = True
    while correct:
        correct = False
        loop_num += 1

        for i in range(1, max_edge_id + 1):
            if i not in edges_cate:
                continue
            cur_cate = edges_cate[i]
            neighbors = edges_adj.get(i, set())

            if cur_cate != 4 and len(neighbors) >= 2:
                num_same = 0
                multi_adj_cate = False
                adj_cate = -1
                num_support = 0

                neighbor_list = sorted(neighbors)
                j = 0
                while j < len(neighbor_list) and not multi_adj_cate:
                    nid = neighbor_list[j]
                    n_cate = edges_cate.get(nid, 4)

                    if n_cate == cur_cate:
                        num_same += 1
                    elif n_cate != 4:
                        adj_cate = n_cate
                        num_support = 1
                        j += 1

                        while not multi_adj_cate and j < len(neighbor_list):
                            nid2 = neighbor_list[j]
                            n2_cate = edges_cate.get(nid2, 4)
                            if n2_cate != 4:
                                if adj_cate != n2_cate:
                                    multi_adj_cate = True
                                else:
                                    num_support += 1
                            j += 1

                        if not multi_adj_cate and num_same == 0 and num_support >= 2:
                            edges_cate[i] = adj_cate
                            correct = True
                            correct_edges += 1
                        break  # exit outer loop after processing
                    j += 1

    return loop_num, correct_edges


def propagation(edges_adj, edges_cate, max_edge_id):
    """Propagation step: propagate classified neighbor categories to uncertain edges.
    Same logic as 3CAC Java Classify3CAC.propagation()
    """
    propagate_edges = 0
    loop_num = 0

    propagate = True
    while propagate:
        propagate = False
        loop_num += 1

        for i in range(1, max_edge_id + 1):
            if i not in edges_cate or edges_cate[i] != 4:
                continue
            neighbors = edges_adj.get(i, set())

            if len(neighbors) >= 1:
                adj_cate = -1
                multi_adj_cate = False

                for nid in neighbors:
                    n_cate = edges_cate.get(nid, 4)
                    if n_cate != 4:
                        if adj_cate == -1:
                            adj_cate = n_cate
                        elif adj_cate != n_cate:
                            multi_adj_cate = True

                if not multi_adj_cate and adj_cate != -1:
                    edges_cate[i] = adj_cate
                    propagate = True
                    propagate_edges += 1

    return loop_num, propagate_edges


def final_contig_category(initial_class, edges_length, edges_cate):
    """Determine final contig category from edge categories.
    For Unicycler: each contig maps to exactly one edge.
    """
    cate_map = {1: 'phage', 2: 'plasmid', 3: 'chromosome', 4: 'uncertain'}
    results = []

    for name, cid, _ in initial_class:
        edge_id = cid  # contig N -> edge N
        if edge_id in edges_cate:
            final_cate = edges_cate[edge_id]
        else:
            final_cate = 4
        results.append((name, cate_map.get(final_cate, 'uncertain')))

    return results


def main():
    parser = argparse.ArgumentParser(description='3CAC-like classifier for Unicycler')
    parser.add_argument('--PPRMeta', required=True, help='PPR-Meta classification CSV')
    parser.add_argument('--plasClass', required=True, help='PlasClass classification output')
    parser.add_argument('--deepVirFinder', required=True, help='DeepVirFinder classification output')
    parser.add_argument('--graph', required=True, help='Assembly graph GFA file')
    parser.add_argument('--output', required=True, help='Output directory')
    args = parser.parse_args()

    # Ensure output directory exists
    outdir = args.output
    if not outdir.endswith('/'):
        outdir += '/'
    os.makedirs(outdir, exist_ok=True)

    # 1. Read PlasClass classification
    plas = read_plasclass(args.plasClass)
    print(f"{len(plas)}  contigs are predicted by plasClass !")

    # 2. Read DeepVirFinder classification
    dvf = read_deepvirfinder(args.deepVirFinder)
    print(f"{len(dvf)}  contigs are predicted by deepVirFinder !")

    # 3. Read PPR-Meta classification
    entries = read_pprmeta(args.PPRMeta)

    # 4. Initial classification
    initial_class, num_phageC, num_plasC, num_uncertain = initial_classification(entries, plas, dvf)

    # Write initial classification
    init_file = outdir + "3CAC_initial_classification.fasta"
    with open(init_file, 'w') as f:
        for name, cid, cate in initial_class:
            f.write(f"{name},{cate}\n")
    print("Finish generating the initial classification of 3CAC !")
    print(f"{num_phageC} contigs classified as phages by PPRMeta were reclassified as chromosome contigs by deepVirFinder.")
    print(f"{num_plasC} contigs classified as plasmids by PPRMeta were reclassified as chromosome contigs by plasClass.")
    print(f"{num_uncertain} contigs classified as phages or plasmids by PPRMeta were reclassified as uncertain.")

    # 5. Read assembly graph
    edges_length, edges_adj, max_edge_id = read_assembly_graph(args.graph)
    print(f"Reading assembly graph {args.graph} !")
    print(f"End Reading assembly graph ! (max_edge_id={max_edge_id})")

    # 6. Initialize edge categories (all uncertain)
    edges_cate = {}
    for eid in range(1, max_edge_id + 1):
        edges_cate[eid] = 4

    # 7. Assign contig categories to edges
    contig_category_to_edges(initial_class, edges_cate)
    print("Edges contained in contigs are assigned with contig categories !")

    # 8. Correction step
    loop_num, correct_edges = correction(edges_adj, edges_cate, max_edge_id)
    print(f"The correction step was repeated for {loop_num} times, and {correct_edges}  classified edges were corrected !")

    # 9. Propagation step
    loop_num, propagate_edges = propagation(edges_adj, edges_cate, max_edge_id)
    print(f"The propagation step was repeated for {loop_num} times, and {propagate_edges}  uncertain edges were classified !")

    # 10. Final contig classification
    results = final_contig_category(initial_class, edges_length, edges_cate)

    # Write final classification
    final_file = outdir + "3CAC_classification.fasta"
    with open(final_file, 'w') as f:
        for name, category in results:
            f.write(f"{name},{category}\n")

    print(f"The classification result of 3CAC can be found in: {final_file}")
    print("Finish 3CAC algorithm !")


if __name__ == '__main__':
    main()
