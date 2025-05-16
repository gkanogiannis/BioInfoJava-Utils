# BioInfoJava-Utils

**BioInfoJava-Utils**  is a modular Java library providing high-performance implementations of core bioinformatics algorithms, such as distance matrix computation and phylogenetic tree construction from **VCF** and **FASTA** files.

This library serves as the computational backend for the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) software suite, which offers a flexible and user-friendly interface to these tools across multiple platforms and environments.

## Integration and Accessibility

The functionality of **BioInfoJava-Utils** is exposed through the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) interface, which is accessible in the following ways:

* ✅ **Python CLI**: through a lightweight [Python wrapper](https://github.com/gkanogiannis/fastreeR/blob/devel/fastreeR.py) that calls the Java backend via `subprocess`
* ✅ **Bioconda**: install with `conda install -c bioconda fastreer`
* ✅ **PyPI**: install with `pip install fastreer`
* ✅ **Docker**: available on [DockerHub](https://hub.docker.com/r/gkanogiannis/fastreer) and [GHCR](https://ghcr.io/gkanogiannis/fastreer) for containerized execution
* ✅ **R / Bioconductor**: via the [`fastreeR`](https://bioconductor.org/packages/fastreeR) package using `rJava`
* ✅ **Pure Java API**: developers can integrate this library directly in Java-based pipelines or software.

## Overview

**BioInfoJava-Utils** provides efficient, scalable, and parallel implementations of widely used bioinformatics algorithms. It is designed for processing large-scale genomic datasets efficiently, supporting both research and production environments.

## Features

* ⚙️ Compute sample-wise **distance matrices** from VCF (cosine) or FASTA (D2S) files
* 🌳 Build **phylogenetic trees** using neighbor-joining algorithm
* 🧬 Support for **hierarchical clustering** with dynamic tree pruning
* 🔄 **Multithreaded** processing for large input files
* 📦 Integrates seamlessly into diverse environments (R, Python, Docker, Java)

## Installation

### Prerequisites

  - Java 8 or higher
  - Maven (for building the project)

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
   mvn clean initialize package
```

This will generate a JAR files in the `bin` directory.

## Usage

The main class for executing the utilities is:

```java
ciat.agrobio.javautils.JavaUtils
```

You can run the utilities via the command line or integrate them into other Java applications.

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
