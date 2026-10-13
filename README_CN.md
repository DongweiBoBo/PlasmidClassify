# Scripts 目录说明文档

本目录包含质粒组装与分类流程（Plasmid Assembly & Classification Pipeline）所需的全部脚本，涵盖 3CAC 分类、质粒序列提取、辅助文件生成以及完整的迭代组装流水线。

---

## 脚本概览

| 脚本名称 | 功能简述 |
|---------|---------|
| `plasmid_circle_pipeline.py` | 完整的迭代质粒组装与分类流水线（主流程） |
| `classify_3CAC_unicycler.py` | 针对 Unicycler 组装结果的 3CAC 分类器 |
| `extract_plasmids_flye.py` | 从 Flye 组装结果中提取质粒序列（区分环状/线性） |
| `extract_plasmids_spades.py` | 从 SPAdes 组装结果中提取质粒序列 |
| `generate_paths.py` | 从 SPAdes GFA 文件生成 3CAC 所需的 `scaffolds.paths` 文件 |

---

## 1. `plasmid_circle_pipeline.py` — 迭代质粒组装与分类流水线

### 功能描述

完整的迭代质粒组装流水线，执行以下流程：

1. 对输入 contigs 进行 BLAST 搜索，比对质粒数据库
2. 根据分类学信息将序列按**属（genus）**分组
3. 对每个属组，使用 minimap2 提取对应的 Nanopore 长读段
4. 使用 **Flye** 进行组装
5. 使用 **PPR_Meta**、**PlasClass**、**DeepVirFinder** 和 **3CAC** 对组装结果进行分类
6. 提取质粒序列，并根据 `assembly_info.txt` 区分环状与非环状质粒
7. 将非环状质粒作为下一轮迭代的输入，重复上述过程

最终输出所有迭代中获得的环状质粒序列。

### 依赖工具

- `blastn`（BLAST+）
- `minimap2`
- `flye`
- `PPR_Meta`
- `classify_fasta.py`（PlasClass）
- `dvf.py`（DeepVirFinder）
- `seqkit`

### 用法

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

### 参数说明

| 参数 | 必填 | 说明 |
|------|------|------|
| `--fasta` | 是 | 输入的 contigs FASTA 文件 |
| `--reads` | 是 | Nanopore 长读段数据（fastq.gz 格式） |
| `--db_index` | 是 | 质粒数据库 BLAST 索引路径 |
| `--tax_dir` | 是 | NCBI 分类学目录（含 `nodes.dmp` 和 `names.dmp`） |
| `--output_dir` | 是 | 输出目录 |
| `--iterations` | 否 | 迭代次数，默认 3 |
| `--threads` | 否 | 线程数，默认 8 |

### 输出文件

- `iteration_N/` — 每次迭代的中间结果目录
- `iteration_N/all_circular_plasmids.fasta` — 当次迭代合并的环状质粒
- `final_circular_plasmids.fasta` — 所有迭代汇总的最终环状质粒结果

---

## 2. `classify_3CAC_unicycler.py` — Unicycler 组装结果的 3CAC 分类

### 功能描述

实现 3CAC（3-source Contig Classification）算法，专门适配 **Unicycler** 组装结果。整合三种分类工具的结果进行综合判定：

- **PPR_Meta**：判断 contig 属于 phage / plasmid / chromosome
- **PlasClass**：对 PPR_Meta 判定为 plasmid 的 contig 进行验证
- **DeepVirFinder**：对 PPR_Meta 判定为 phage 的 contig 进行验证

分类完成后，利用组装图（GFA）进行两步图矫正：
- **Correction**：根据邻居一致性修正错误分类的 edge
- **Propagation**：将已分类 edge 的类别传播到不确定（uncertain）的 edge

最终将 edge 分类结果映射回 contig，输出最终分类。

### 分类类别

| 编号 | 类别 |
|------|------|
| 1 | phage |
| 2 | plasmid |
| 3 | chromosome |
| 4 | uncertain |

### 用法

```bash
python3 classify_3CAC_unicycler.py \
    --PPRMeta Sample/PPR_Meta.csv \
    --plasClass Sample/PlasClass.probs.out \
    --deepVirFinder Sample/dvfpred.txt \
    --graph Sample/assembly.gfa \
    --output Sample/3CAC_output/
```

### 参数说明

| 参数 | 必填 | 说明 |
|------|------|------|
| `--PPRMeta` | 是 | PPR-Meta 分类结果 CSV 文件 |
| `--plasClass` | 是 | PlasClass 分类输出文件 |
| `--deepVirFinder` | 是 | DeepVirFinder 分类输出文件 |
| `--graph` | 是 | 组装图 GFA 文件 |
| `--output` | 是 | 输出目录 |

### 输出文件

- `3CAC_initial_classification.fasta` — 初始分类结果（未经图矫正）
- `3CAC_classification.fasta` — 最终分类结果（经 Correction + Propagation）

每行格式为 `contig_name,category`，其中 category 为 `phage`/`plasmid`/`chromosome`/`uncertain`。

---

## 3. `extract_plasmids_flye.py` — 从 Flye 组装结果提取质粒序列

### 功能描述

读取 3CAC 分类结果，识别被分类为 **plasmid** 的 contig，再根据 Flye 的 `assembly_info.txt` 中的 **circ.**（circularity）列区分环状与非环状质粒，分别输出到两个 FASTA 文件。

### 用法

```bash
python3 extract_plasmids_flye.py \
    --classification Sample.flye/3CAC_classification.fasta \
    --assembly_info Sample.flye/assembly_info.txt \
    --assembly_fasta Sample.flye/assembly.fasta \
    --output_dir Sample.flye/
```

### 参数说明

| 参数 | 必填 | 说明 |
|------|------|------|
| `--classification` | 是 | 3CAC 分类结果文件 |
| `--assembly_info` | 是 | Flye 的 `assembly_info.txt` |
| `--assembly_fasta` | 是 | Flye 原始组装 FASTA 文件 |
| `--output_dir` | 是 | 输出目录 |

### 输出文件

- `plasmid_circular.fasta` — 环状质粒序列
- `plasmid_linear.fasta` — 非环状质粒序列

---

## 4. `extract_plasmids_spades.py` — 从 SPAdes 组装结果提取质粒序列

### 功能描述

读取 3CAC 分类结果，从 SPAdes 组装的 FASTA 文件中提取被分类为 **plasmid** 的 contig 序列。通过解析 FASTA header 中的 `NODE_X` 标识与分类结果进行匹配。

> **注意**：SPAdes 组装结果不区分环状/线性，因此仅输出一个质粒 FASTA 文件。

### 用法

```bash
python3 extract_plasmids_spades.py \
    --classification Sample.spades/3CAC_classification.fasta \
    --assembly_fasta Sample.spades/scaffolds.fasta \
    --output Sample.spades/plasmid.fasta
```

### 参数说明

| 参数 | 必填 | 说明 |
|------|------|------|
| `--classification` | 是 | 3CAC 分类结果文件 |
| `--assembly_fasta` | 是 | SPAdes 组装 FASTA 文件 |
| `--output` | 是 | 输出 FASTA 文件路径 |

### 输出文件

- 指定的输出 FASTA 文件，包含所有被分类为 plasmid 的 contig 序列

---

## 5. `generate_paths.py` — 生成 3CAC 所需的 scaffolds.paths 文件

### 功能描述

3CAC 算法在处理 SPAdes 组装结果时需要 `scaffolds.paths` 文件来描述 contig 与 edge 的对应关系。此脚本从 SPAdes 输出的 GFA 文件和 FASTA 文件中自动提取 edge ID 与 contig header 的映射关系，生成符合 3CAC 格式要求的 `scaffolds.paths` 文件。

每条 contig 生成 4 行记录：

```
contig_header
edge_id+
contig_header'
edge_id-
```

### 用法

```bash
python3 generate_paths.py <gfa_file> <fasta_file> <output_file>
```

**示例：**

```bash
python3 generate_paths.py \
    SPAdes/assembly_graph_with_scaffolds.gfa \
    SPAdes/scaffolds.fasta \
    SPAdes/scaffolds.paths
```

### 参数说明

| 参数 | 位置 | 说明 |
|------|------|------|
| `gfa_file` | 1 | SPAdes 输出的 GFA 文件 |
| `fasta_file` | 2 | SPAdes 输出的 scaffolds FASTA 文件 |
| `output_file` | 3 | 输出的 `scaffolds.paths` 文件路径 |

---

## 典型工作流程

以下是一个样本从组装到质粒提取的完整流程示例：

### Flye 组装流程

```bash
# 1. Flye 组装（略）

# 2. 运行 PPR_Meta, PlasClass, DeepVirFinder（略）

# 3. 3CAC 分类（使用 Unicycler 版分类器，也适用于 Flye GFA）
python3 scripts/classify_3CAC_unicycler.py \
    --PPRMeta SYZ.0.flye/SYZ.0.PPR_Meta.csv \
    --plasClass SYZ.0.flye/SYZ.0.PlasClass.probs.out \
    --deepVirFinder SYZ.0.flye/dvfpred.txt \
    --graph SYZ.0.flye/assembly_graph.gfa \
    --output SYZ.0.flye/

# 4. 提取质粒序列
python3 scripts/extract_plasmids_flye.py \
    --classification SYZ.0.flye/3CAC_classification.fasta \
    --assembly_info SYZ.0.flye/assembly_info.txt \
    --assembly_fasta SYZ.0.flye/assembly.fasta \
    --output_dir SYZ.0.flye/
```

### SPAdes 组装流程

```bash
# 1. SPAdes 组装（略）

# 2. 生成 scaffolds.paths
python3 scripts/generate_paths.py \
    SPAdes/assembly_graph_with_scaffolds.gfa \
    SPAdes/scaffolds.fasta \
    SPAdes/scaffolds.paths

# 3. 运行 PPR_Meta, PlasClass, DeepVirFinder（略）

# 4. 3CAC 分类（略）

# 5. 提取质粒序列
python3 scripts/extract_plasmids_spades.py \
    --classification SYZ.0.spades/3CAC_classification.fasta \
    --assembly_fasta SYZ.0.spades/scaffolds.fasta \
    --output SYZ.0.spades/plasmid.fasta
```

### 完整迭代流水线

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
