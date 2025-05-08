# BioInfoJava-Utils

**BioInfoJava-Utils** is a collection of Java-based utilities designed to facilitate various bioinformatics analyses, including distance calculations, phylogenetic tree construction, and hierarchical clustering on VCF and FASTA files.

## Overview

This toolkit provides efficient and parallelized implementations of common bioinformatics algorithms, making it suitable for large-scale genomic data analyses. The utilities are also integrated into the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) R [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/fastreeR.html) package via `rJava` and the relevant [`fastreeR.py`](https://github.com/gkanogiannis/fastreeR/bin) CLI wrapper, allowing seamless use within R, Linux and Windows console environments.

## Features

- Compute distance matrices from VCF and FASTA files.
- Construct phylogenetic trees using agglomerative neighbor-joining methods.
- Perform hierarchical clustering and dynamic tree pruning.
- Parallel processing capabilities for handling large datasets.

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
   mvn clean install package
   ```

This will generate a JAR file in the `target` directory.

## Usage

The main class for executing the utilities is:

```java
ciat.agrobio.javautils.JavaUtils
```

You can run the utilities via the command line or integrate them into other Java applications.

## Integration with R

The functionalities provided by BioInfoJava-Utils are accessible in R through the [`fastreeR`](https://github.com/gkanogiannis/fastreeR) package, which leverages `rJava` for interfacing with Java code. This integration allows users to perform complex bioinformatics analyses directly within R.

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
