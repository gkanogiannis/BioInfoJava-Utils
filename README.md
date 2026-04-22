# BioInfoJava-Utils

**BioInfoJava-Utils** is a modular Java library providing high-performance implementations of core bioinformatics algorithms, such as distance matrix computation and phylogenetic tree construction (via hierarchical clustering) from **VCF** and **FASTA** files.

This library serves as the computational backend for the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) software suite, which offers a flexible and user-friendly interface to these tools across multiple platforms and environments.

## Integration and Accessibility

The functionality of **BioInfoJava-Utils** is exposed through the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) interface, which is accessible in the following ways:

* 🆕 **Java Backend ([v2.7.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.7.0)) !!** introduces **windowed / streaming VCF distance & tree output**. Emit one distance matrix (or Newick tree) per genomic window of N base pairs (`--window-bp`) or per N consecutive variants (`--window-variants`) for `VCF2DIST` and `VCF2TREE`, with optional long-form TSV output (`--long`).
* Java Backend ([v2.5.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.5.0)) introduces **embedding-based distance calculation** for VCF files. Provide pre-computed variant embeddings (from genomic language models like [BioFM](https://huggingface.co/m42-health/BioFM-265M), DNA-BERT, Nucleotide Transformer, etc.) to weight variant contributions during distance computation.
* Java Backend ([v2.3.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.3.0)) supports reading from gzip (for example .gz), bzip2 (for example .bz2) and xz compressed VCF files.
* Java Backend ([v2.2.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.2.0)) implements streaming bootstrap; from VCF file get a newick tree with encoded bootstrap support values
* Java Backend ([v2.0.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/2.0.0)) 100x times **FAST**re**ER** and only a couple hundred MB RAM needed. Java 11+ suggested.
* **Bioconda**: install with `conda install -c bioconda fastreer` ([recipe](https://bioconda.github.io/recipes/fastreer/README.html))
* **Docker**: available on [DockerHub](https://hub.docker.com/r/gkanogiannis/fastreer) and [GHCR](https://ghcr.io/gkanogiannis/fastreer) for containerized execution
* **PyPI**: install with `pip install fastreer` ([repository](https://pypi.org/project/fastreer/))
* **Python CLI**: through a lightweight [Python wrapper](https://github.com/gkanogiannis/fastreeR/blob/devel/fastreeR.py) that calls the Java backend
* **R / Bioconductor**: via `rJava` ([package](https://bioconductor.org/packages/fastreeR/))
* **Galaxy**: available on Galaxy [Toolshed](https://toolshed.g2.bx.psu.edu/view/gkanogiannis/fastreer/26013530719e).
* **Pure Java API**: developers can integrate this library directly in Java-based pipelines or software.

## Overview

**BioInfoJava-Utils** provides efficient, scalable, and parallel implementations of widely used bioinformatics algorithms. It is designed for processing large-scale genomic datasets efficiently, supporting both research and production environments.

## Features

* Reads directly from plain, gzip, bzip2 or xz VCF files.
* 🪟 **Windowed / streaming output** emits one distance matrix or Newick tree per genomic window (by base pairs or variant count) for `VCF2DIST` and `VCF2TREE`
* 🧠 **Embedding-based distance calculation** using pre-computed variant embeddings from genomic language models
* 🥾 **Streaming bootstrap** support in the VCF2TREE utility
* 🚀 Ultra-fast with a superior multithreaded concurrency model
  and minimal RAM usage, **from GBs down to just MBs!**
* ⚙️ Compute sample-wise **distance matrices** from VCF (cosine) or FASTA (D2S) files
* 🌳 Build **phylogenetic trees** using **hierarchical clustering** (single, complete, or average linkage; complete by default)
* 🧬 Support for **hierarchical clustering** with dynamic tree pruning
* 🔄 **Multithreaded** processing for large input files
* 📦 Integrates seamlessly into diverse environments (R, Python, Docker, Java)

## Installation

### Prerequisites

* Java 11 or higher
* Maven (for building the project)

### Building from Source

1. Clone the repository:

```bash
   git clone https://github.com/gkanogiannis/BioInfoJava-Utils.git
```

2. Navigate to the project directory:

```bash
   cd BioInfoJava-Utils
```

3. Build the project using Maven:

```bash
   mvn clean package install
```

This will generate a JAR files in the `bin` directory.

## Usage

The main class for executing the utilities is:

```java
com.gkano.bioinfo.javautils.JavaUtils
```

You can run the utilities via the command line or integrate them into other Java applications.

```bash
java -jar bin/BioInfoJavaUtils-VERSION-jar-with-dependencies.jar --help
```

## Embedding-Based Distance Calculation

Version 2.5.0 introduces support for **embedding-based distance calculation** in `VCF2DIST` and `VCF2TREE`. This feature allows you to incorporate pre-computed variant embeddings (e.g., from genomic language models like [BioFM](https://huggingface.co/m42-health/BioFM-265M), DNA-BERT, Nucleotide Transformer, or custom embeddings) to compute distances in embedding space rather than genotype space.

### How It Works

Instead of computing cosine similarity directly from genotype vectors, the embedding mode:

1. Projects each sample into embedding space: `H_i = Σ_v dosage_i^v × e_v`
2. Computes cosine distance between sample embeddings

This captures functional relationships between variants - samples with alleles at functionally similar positions become more similar in embedding space.

### Embedding File Formats

**TSV Format:**

```tsv
#VARIANT_ID  DIM_0   DIM_1   DIM_2   ...
chr1:12345:A:G  0.123   -0.456  0.789   ...
chr1:67890:C:T  0.567   0.123   -0.890  ...
```

**HuggingFace JSON Format:**

```json
{
  "model_name": "genomic-model-name",
  "embedding_dim": 768,
  "variants": [
    {"id": "chr1:12345:A:G", "embedding": [0.123, -0.456, ...]},
    {"id": "chr1:67890:C:T", "embedding": [0.567, 0.123, ...]}
  ]
}
```

### Command Line Options

| Option                 | Description                                                                     |
|------------------------|---------------------------------------------------------------------------------|
| `-e, --embeddings`     | Path to variant embeddings file                                                 |
| `--embeddings-format`  | Format: `TSV` or `HUGGINGFACE` (auto-detected if not specified)                 |
| `--variant-key`        | Variant key format: `CHROM_POS`, `CHROM_POS_REF_ALT` (default), or `VCF_ID`     |

### Examples

```bash
# Distance matrix with embeddings (TSV format, auto-detected)
java -jar BioInfoJavaUtils.jar VCF2DIST \
  -i samples.vcf.gz \
  -o distances.tsv \
  -e variant_embeddings.tsv \
  -t 4

# Tree with embeddings and bootstrap (HuggingFace format)
java -jar BioInfoJavaUtils.jar VCF2TREE \
  -i samples.vcf.gz \
  -o tree.nwk \
  -e embeddings.json \
  --embeddings-format HUGGINGFACE \
  -b 100

# Standard mode (no embeddings) - existing behavior
java -jar BioInfoJavaUtils.jar VCF2DIST \
  -i samples.vcf.gz \
  -o distances.tsv
```

Variants without matching embeddings are automatically skipped, and the tool reports how many variants were used vs. skipped.

## Windowed / Streaming Output

Version 2.7.0 introduces **windowed output** for `VCF2DIST` and `VCF2TREE`. Instead of producing a single genome-wide distance matrix or tree, the tools can stream one matrix (or Newick tree) per genomic window. This enables local-ancestry analyses, introgression scans, recombination-rate studies, and any workflow that needs sample relationships measured along the genome.

### How It Works

Variants are streamed in input order and grouped into windows defined either by base-pair span (`--window-bp`) or by consecutive variant count (`--window-variants`). When a window closes, all worker threads synchronize on a barrier, the per-window distance matrix is reduced from shared accumulators, the writer emits the window, and the accumulators are zeroed before the next window opens. Windows never straddle chromosomes; a contig change always closes the current window.

The non-windowed code path is unchanged and remains byte-identical to previous releases.

### Command Line Options

| Option              | Description                                                                                                |
|---------------------|------------------------------------------------------------------------------------------------------------|
| `--window-bp`       | Emit one matrix/tree per window of N base pairs (mutually exclusive with `--window-variants`)              |
| `--window-variants` | Emit one matrix/tree per N consecutive variants (mutually exclusive with `--window-bp`)                    |
| `--step`            | Window step. Defaults to window size (tiled). Sliding windows (`step != size`) are not yet implemented.    |
| `--min-variants`    | Minimum number of variants required to emit a window (default 1; smaller windows are skipped silently)     |
| `--long`            | (`VCF2DIST` only) Emit long-form TSV `chrom, start, end, sample_i, sample_j, dist` instead of matrices     |

### Output Formats

`VCF2DIST` default (concatenated matrices); one block per window:

```text
# window chrom=chr1 start=0 end=100000 nvariants=842 nsamples=3
3	842
s1	0	0.4231	0.5102
s2	0.4231	0	0.3987
s3	0.5102	0.3987	0
# window chrom=chr1 start=100000 end=200000 nvariants=917 nsamples=3
...
```

`VCF2DIST --long`; single TSV with one row per sample pair per window:

```text
chrom	start	end	sample_i	sample_j	dist
chr1	0	100000	s1	s2	0.4231
chr1	0	100000	s1	s3	0.5102
chr1	0	100000	s2	s3	0.3987
...
```

`VCF2TREE`; one Newick tree per window, prefixed by a header comment:

```text
# window chrom=chr1 start=0 end=100000 nvariants=842 nsamples=3
(s1:0.21,(s2:0.19,s3:0.18):0.05);
# window chrom=chr1 start=100000 end=200000 nvariants=917 nsamples=3
(s2:0.20,(s1:0.22,s3:0.17):0.04);
...
```

### Examples

```bash
# Distance matrices in 100kb tiled windows
java -jar BioInfoJavaUtils.jar VCF2DIST \
  -i samples.vcf.gz \
  -o per_window.dist \
  --window-bp 100000 \
  -t 4

# Long-form TSV, one matrix per 500 consecutive variants
java -jar BioInfoJavaUtils.jar VCF2DIST \
  -i samples.vcf.gz \
  -o per_window.tsv \
  --window-variants 500 \
  --long \
  -t 4

# Per-window phylogenetic trees (Newick)
java -jar BioInfoJavaUtils.jar VCF2TREE \
  -i samples.vcf.gz \
  -o per_window.nwk \
  --window-bp 250000 \
  -t 4

# Skip windows with fewer than 50 variants
java -jar BioInfoJavaUtils.jar VCF2DIST \
  -i samples.vcf.gz \
  -o per_window.dist \
  --window-bp 100000 \
  --min-variants 50
```

### Limitations

- **Sliding windows** (`--step` different from window size) are reserved for a future release; passing them throws an error.
- **Bootstrap** (`-b` / `--bootstrap`) is rejected when combined with windowing.
- **Embeddings** (`-e` / `--embeddings`) are rejected when combined with windowing.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

## Citation

If you use BioInfoJava-Utils in your research, please cite the following:

```text
Gkanogiannis, A. et al. A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes. BMC Bioinformatics 17, 311 (2016). https://doi.org/10.1186/s12859-016-1186-3
```

## Author

**Anestis Gkanogiannis**  
Bioinformatics/ML Scientist  
Linkedin: [https://www.linkedin.com/in/anestis-gkanogiannis/](https://www.linkedin.com/in/anestis-gkanogiannis/)  
Website: [https://www.gkanogiannis.com](https://www.gkanogiannis.com)  
ORCID: [0000-0002-6441-0688](https://orcid.org/0000-0002-6441-0688)
