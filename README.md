# BioInfoJava-Utils

**BioInfoJava-Utils** is a modular Java library providing high-performance implementations of core bioinformatics algorithms, such as distance matrix computation and phylogenetic tree construction from **VCF** and **FASTA** files.

This library serves as the computational backend for the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) software suite, which offers a flexible and user-friendly interface to these tools across multiple platforms and environments.

## Integration and Accessibility

The functionality of **BioInfoJava-Utils** is exposed through the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) interface, which is accessible in the following ways:

* 🆕 **Java Backend ([v2.2.0](https://github.com/gkanogiannis/BioInfoJava-Utils/releases/tag/v2.2.0)) !!** implements streaming bootstrap; from VCF file get a newick tree with encoded bootstrap support values
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

```java
java -jar bin/BioInfoJavaUtils-VERSION-jar-with-dependencies.jar --help
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

## Citation

If you use BioInfoJava-Utils in your research, please cite the following:

```
Gkanogiannis, A. et al. A scalable assembly-free variable selection algorithm for biomarker discovery from metagenomes. BMC Bioinformatics 17, 311 (2016). https://doi.org/10.1186/s12859-016-1186-3
```

## Author

**Anestis Gkanogiannis**  
Bioinformatics/ML Scientist  
Website: [https://www.gkanogiannis.com](https://www.gkanogiannis.com)  
ORCID: [0000-0002-6441-0688](https://orcid.org/0000-0002-6441-0688)
