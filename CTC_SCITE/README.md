# CTC-SCITE
## Description

**CTC-SCITE** is a software package to infer cell lineages from single-cells and clusters of single-cells.
It can be used to infer the clonality of circulating tumor cell (CTC) clusters.
CTCs are often shed into the bloodstream as clusters of two or more physically associated cells and frequently admix
with tumor-associated cells like white blood cells. When these clusters cannot be physically
broken into single cells without harming cell integrity, the joint DNA material of all cells in the cluster
is sequenced in a single sample, resulting in an aggregate set of variants of the CTC cluster.
Given these aggregate SNV read count profiles of clusters of single cells, the algorithm jointly infers the genotypes of
the cells and their genealogical relationships.
By splitting the constituent cells in the tree, it can infer each cell's individual genotype
allowing us to the the mixed and assess clonality of the cell clusters.


## Availability

**CTC-SCITE** is freely available under a GPL3 license at http://CTC-SCITE


##    Installation

**CTC-SCITE** is available as a conda package from Bioconda for linux-64, linux-aarch64, osx-64 and osx-arm64 architectures. To install **CTC-SCITE**, please make sure you have conda installed on your machine, e.g. following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Once conda is installed, you can install **CTC-SCITE** into an environment `env_name` by running the command

```bash
conda create -n env_name -c bioconda ctc-scite
```
in the terminal.

To run the code, activate the conda environment:

```bash
conda activate env_name
```

To check of the installation succeeded run

```bash
CTC-SCITE -help 
```
This should print a help message explaining the command line interface.

In case you are a user of Windows 10 or higher, you can run **CTC-SCITE** by installing wsl (Windows subsystem for Linux). Open a Powershell terminal and enter

```Powershell
wsl --install
```
Please note that the installation of wsl may take a while. Now launch the Linux subsystem with

```Powershell
wsl
```
and follow the installation intructions for Linux users as above.

##  Input Files

### 1. Mutation Matrix

The mutation input data consists of a matrix in tab-separated format. Each row corresponds to an SNV site, and pair of columns to a cell or cell cluster (sample).
The entries for each mutation row and each sample column pair are reference read counts in the first of the columns, followed by alternative read counts in the second of the columns.
Each column specifies the mutation profile of a single cell, and each row
represents one mutation. Columns are separated by a white space character. The first columns of the file specify mutations information (genomic position, reference and alternative allele).

An example file is given as input/ex.txt.

### 2. A sample description file.

A .tsv file specifying the metadata for each of the samples. Each row corresponds to a sample. The first column specifies the sample name, the second column contains the total number of cells in the sample,
the third column specifies the total number of tumor cells in the sample and the fourth column specifies the number of White blood cells in the sample.

An exaple file is given as input/ex_samples_nodeDescription.tsv

## Run CTC-SCTE

**CTC-SCITE** can then be executed as follows:

```bash
CTC-SCITE -i input/ex.txt -samples input/ex_samples_nodeDescription.tsv -r 1 -l 1000 -g 1 -e 0.1 -p 1 -o input/test
```

##  Output Files

### 1. ML/MAP trees

ML/MAP trees are written to files in GraphViz format. Files are numbered consecutively (e.g. `ex_rank_1.gv`, `ex_rank_2.gv`, ...). The base name of the output file is derived from the name of the input file (unless a different name is specified via `-o <filename>`).

### 2. Samples from the posterior distribution

When the `-p <INT>` option is set, **CTC-SCITE** samples from the posterior distribution, and writes the sampled trees (i.e. as Prüfer sequence) together with their scores and learned error rates to a single tab-separated file (one sample per line). The name of the output file is derived from the input file name using the ending `_post_sampling.tsv`.

To enable sampling from the posterior distribution, set `-p <INT>`. **CTC-SCITE**, will write the sampled trees to a file using the parent vector format. The value of <INT> specifies how dense the sampling is. The name of the output file is derived from the input file name using the ending `.postSampling.tsv`. To make sure that SCITE samples from the posterior distribution -p <INT> needs to be combined with `-g 1` (gamma is set to 1).





## Parameters

*	`-i <filename>`  Replace \<filename\> with the file containing the mutation matrix

*	`-samples <sampledesciption>`  Replace \<sampledescription\> by the file containing the sample annotations

* 	`-r <INT>`  Set \<INT\> to the desired number of repetitions of the MCMC

* 	`-l <INT>`  Set \<INT\> to the desired chain length of each MCMC repetition

*   `-g <DOUBLE>`  For ML/MAP computation only: Set \<DOUBLE\> to the desired value of gamma (gamma > 1: more local exploration, possibly local optimum; gamma < 1: easier to explore the space, but less deeply). The default value of gamma is 1. This is necessary for the MCMC chain to guarantee asymptotic convergence to the posterior distribution.

*   `-e <DOUBLE>`  Invokes the learning of error rates. Set \<DOUBLE\> to a value between zero and one to specify the probability to chose the move for changing the error rate in the MCMC

*   `-p <INT>`  When setting this option, CTC-SCITE samples from the posterior distribution, and writes the trees to a file using the parent vector format. The value of <INT> specifies how dense the sampling is. In this case the parameter -g must be set to 1

*   `-o <filename>` Optional. Replace <filename> with the desired base of the output file to overwrite the default output file names
