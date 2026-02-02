# BioInfoJava-Utils

**BioInfoJava-Utils** is a modular Java library providing high-performance implementations of core bioinformatics algorithms, such as distance matrix computation and phylogenetic tree construction from **VCF** and **FASTA** files.

This library serves as the computational backend for the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) software suite, which offers a flexible and user-friendly interface to these tools across multiple platforms and environments.

## Integration and Accessibility

The functionality of **BioInfoJava-Utils** is exposed through the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) interface, which is accessible in the following ways:

* 🆕 **Java Backend ([v2.5.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.5.0)) !!** introduces **embedding-based distance calculation** for VCF files. Provide pre-computed variant embeddings (from genomic language models like [BioFM](https://huggingface.co/m42-health/BioFM-265M), DNA-BERT, Nucleotide Transformer, etc.) to weight variant contributions during distance computation.
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
* 🧠 **Embedding-based distance calculation** using pre-computed variant embeddings from genomic language models
* 🥾 **Streaming bootstrap** support in the VCF2TREE utility
* 🚀 Ultra-fast with a superior multithreaded concurrency model
  and minimal RAM usage, **from GBs down to just MBs!**
* ⚙️ Compute sample-wise **distance matrices** from VCF (cosine) or FASTA (D2S) files
* 🌳 Build **phylogenetic trees** using neighbor-joining algorithm
* 🧬 Support for **hierarchical clustering** with dynamic tree pruning
* 🔄 **Multithreaded** processing for large input files
* 📦 Integrates seamlessly into diverse environments (R, Python, Docker, Java)

## Installation

### Prerequisites

* Java 17 or higher
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
