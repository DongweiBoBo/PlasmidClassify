# Scripts Directory Documentation

This directory contains all scripts for the Plasmid Assembly and Classification Pipeline, covering 3CAC classification, plasmid sequence extraction, auxiliary file generation, and the complete iterative assembly pipeline.

---

## Script Overview

| Script | Description |
|--------|-------------|
| `plasmid_circle_pipeline.py` | Complete iterative plasmid assembly and classification pipeline (main workflow) |
| `classify_3CAC_unicycler.py` | 3CAC classifier adapted for Unicycler assembly results |
| `extract_plasmids_flye.py` | Extract plasmid sequences from Flye assembly results (circular/linear separation) |
| `extract_plasmids_spades.py` | Extract plasmid sequences from SPAdes assembly results |
| `generate_paths.py` | Generate `scaffolds.paths` file required by 3CAC from SPAdes GFA output |

---

## 1. `plasmid_circle_pipeline.py` — Iterative Plasmid Assembly & Classification Pipeline

### Description

A complete iterative plasmid assembly pipeline that performs the following workflow:

1. BLAST search of input contigs against a plasmid database
2. Group sequences by **genus** based on taxonomy information
3. For each genus group, extract corresponding Nanopore long reads using minimap2
4. Assemble reads using **Flye**
5. Classify assemblies using **PPR_Meta**, **PlasClass**, **DeepVirFinder**, and **3CAC**
6. Extract plasmid sequences and separate circular from non-circular plasmids based on `assembly_info.txt`
7. Feed non-circular plasmids into the next iteration

The final output merges all circular plasmids recovered across all iterations.

### Dependencies

- `blastn` (BLAST+)
- `minimap2`
- `flye`
- `PPR_Meta`
- `classify_fasta.py` (PlasClass)
- `dvf.py` (DeepVirFinder)
- `seqkit`

### Usage

```bash
python3 plasmid_circle_pipeline.py \
    --fasta input.fasta \
    --reads reads.fastq.gz \
    --db_index /path/to/plasmid_db.ndb \
    --tax_dir /path/to/taxonomy \
    --output_dir output \
    --iterations 3 \
    --threads 8
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--fasta` | Yes | Input contigs FASTA file |
| `--reads` | Yes | Nanopore long-read data (fastq.gz format) |
| `--db_index` | Yes | Plasmid database BLAST index path |
| `--tax_dir` | Yes | NCBI taxonomy directory (containing `nodes.dmp` and `names.dmp`) |
| `--output_dir` | Yes | Output directory |
| `--iterations` | No | Number of iterations (default: 3) |
| `--threads` | No | Number of threads (default: 8) |

### Output Files

- `iteration_N/` — Intermediate results directory for each iteration
- `iteration_N/all_circular_plasmids.fasta` — Merged circular plasmids from that iteration
- `final_circular_plasmids.fasta` — Final merged circular plasmids from all iterations

---

## 2. `classify_3CAC_unicycler.py` — 3CAC Classification for Unicycler Assemblies

### Description

Implements the 3CAC (3-source Contig Classification) algorithm, specifically adapted for **Unicycler** assembly results. It integrates predictions from three classification tools:

- **PPR_Meta**: Classifies contigs as phage / plasmid / chromosome
- **PlasClass**: Validates contigs classified as plasmid by PPR_Meta
- **DeepVirFinder**: Validates contigs classified as phage by PPR_Meta

After initial classification, two graph-based refinement steps are performed using the assembly graph (GFA):
- **Correction**: Corrects misclassified edges based on neighbor consensus
- **Propagation**: Propagates classified edge categories to uncertain edges

The final edge classifications are mapped back to contigs.

### Classification Categories

| Code | Category |
|------|----------|
| 1 | phage |
| 2 | plasmid |
| 3 | chromosome |
| 4 | uncertain |

### Usage

```bash
python3 classify_3CAC_unicycler.py \
    --PPRMeta Sample/PPR_Meta.csv \
    --plasClass Sample/PlasClass.probs.out \
    --deepVirFinder Sample/dvfpred.txt \
    --graph Sample/assembly.gfa \
    --output Sample/3CAC_output/
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--PPRMeta` | Yes | PPR-Meta classification CSV file |
| `--plasClass` | Yes | PlasClass classification output file |
| `--deepVirFinder` | Yes | DeepVirFinder classification output file |
| `--graph` | Yes | Assembly graph GFA file |
| `--output` | Yes | Output directory |

### Output Files

- `3CAC_initial_classification.fasta` — Initial classification results (before graph correction)
- `3CAC_classification.fasta` — Final classification results (after Correction + Propagation)

Each line follows the format `contig_name,category`, where category is `phage`, `plasmid`, `chromosome`, or `uncertain`.

---

## 3. `extract_plasmids_flye.py` — Extract Plasmid Sequences from Flye Assembly

### Description

Reads 3CAC classification results to identify contigs classified as **plasmid**, then separates them into circular and non-circular plasmids based on the **circ.** (circularity) column in Flye's `assembly_info.txt`. Outputs are written to two separate FASTA files.

### Usage

```bash
python3 extract_plasmids_flye.py \
    --classification Sample.flye/3CAC_classification.fasta \
    --assembly_info Sample.flye/assembly_info.txt \
    --assembly_fasta Sample.flye/assembly.fasta \
    --output_dir Sample.flye/
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--classification` | Yes | 3CAC classification result file |
| `--assembly_info` | Yes | Flye `assembly_info.txt` file |
| `--assembly_fasta` | Yes | Original Flye assembly FASTA file |
| `--output_dir` | Yes | Output directory |

### Output Files

- `plasmid_circular.fasta` — Circular plasmid sequences
- `plasmid_linear.fasta` — Non-circular plasmid sequences

---

## 4. `extract_plasmids_spades.py` — Extract Plasmid Sequences from SPAdes Assembly

### Description

Reads 3CAC classification results and extracts contigs classified as **plasmid** from SPAdes assembly FASTA files. Matches are performed by parsing the `NODE_X` identifier in FASTA headers against the classification results.

> **Note**: SPAdes assembly results do not distinguish circular from linear sequences, so only a single plasmid FASTA file is produced.

### Usage

```bash
python3 extract_plasmids_spades.py \
    --classification Sample.spades/3CAC_classification.fasta \
    --assembly_fasta Sample.spades/scaffolds.fasta \
    --output Sample.spades/plasmid.fasta
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--classification` | Yes | 3CAC classification result file |
| `--assembly_fasta` | Yes | SPAdes assembly FASTA file |
| `--output` | Yes | Output FASTA file path |

### Output Files

- The specified output FASTA file containing all contigs classified as plasmid

---

## 5. `generate_paths.py` — Generate scaffolds.paths File for 3CAC

### Description

The 3CAC algorithm requires a `scaffolds.paths` file when processing SPAdes assembly results to describe the mapping between contigs and edges. This script automatically extracts edge IDs and contig header mappings from SPAdes GFA and FASTA files, generating a `scaffolds.paths` file in the format required by 3CAC.

Each contig generates 4 lines:

```
contig_header
edge_id+
contig_header'
edge_id-
```

### Usage

```bash
python3 generate_paths.py <gfa_file> <fasta_file> <output_file>
```

**Example:**

```bash
python3 generate_paths.py \
    SPAdes/assembly_graph_with_scaffolds.gfa \
    SPAdes/scaffolds.fasta \
    SPAdes/scaffolds.paths
```

### Parameters

| Parameter | Position | Description |
|-----------|----------|-------------|
| `gfa_file` | 1 | SPAdes GFA output file |
| `fasta_file` | 2 | SPAdes scaffolds FASTA file |
| `output_file` | 3 | Output `scaffolds.paths` file path |

---

## Typical Workflows

### Flye Assembly Workflow

```bash
# 1. Flye assembly (omitted)

# 2. Run PPR_Meta, PlasClass, DeepVirFinder (omitted)

# 3. 3CAC classification
python3 scripts/classify_3CAC_unicycler.py \
    --PPRMeta SYZ.0.flye/SYZ.0.PPR_Meta.csv \
    --plasClass SYZ.0.flye/SYZ.0.PlasClass.probs.out \
    --deepVirFinder SYZ.0.flye/dvfpred.txt \
    --graph SYZ.0.flye/assembly_graph.gfa \
    --output SYZ.0.flye/

# 4. Extract plasmid sequences
python3 scripts/extract_plasmids_flye.py \
    --classification SYZ.0.flye/3CAC_classification.fasta \
    --assembly_info SYZ.0.flye/assembly_info.txt \
    --assembly_fasta SYZ.0.flye/assembly.fasta \
    --output_dir SYZ.0.flye/
```

### SPAdes Assembly Workflow

```bash
# 1. SPAdes assembly (omitted)

# 2. Generate scaffolds.paths
python3 scripts/generate_paths.py \
    SPAdes/assembly_graph_with_scaffolds.gfa \
    SPAdes/scaffolds.fasta \
    SPAdes/scaffolds.paths

# 3. Run PPR_Meta, PlasClass, DeepVirFinder (omitted)

# 4. 3CAC classification (omitted)

# 5. Extract plasmid sequences
python3 scripts/extract_plasmids_spades.py \
    --classification SYZ.0.spades/3CAC_classification.fasta \
    --assembly_fasta SYZ.0.spades/scaffolds.fasta \
    --output SYZ.0.spades/plasmid.fasta
```

### Complete Iterative Pipeline

```bash
python3 scripts/plasmid_circle_pipeline.py \
    --fasta combine/SYZ.0.mix.fasta \
    --reads nanopore/SYZ.0.flye.fasta \
    --db_index /path/to/plasmid_db.ndb \
    --tax_dir /path/to/taxonomy \
    --output_dir SYZ.0.pipeline/ \
    --iterations 3 \
    --threads 16
```
