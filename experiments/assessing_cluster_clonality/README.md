[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

# Snakemake workflow for full reproducibility of the project
============================================================

## Description
--------------

This directory contains the code executed to perform the interrogation of clonality of CTC clusters. To facilitate full reproducibility, it has been implemented as a **snakemake** workflow. Some code for one-time manual manipulation of the data can be found in the subdirectory `sandbox/`.


## Usage
--------

The pipeline is implemented in snakemake 8. For the installation, please follow the [official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html 'Snakemake docs').

To rerun the analysis, navigate to this directory in the terminal and run

`snakemake --configfile config/config.yaml --cores 4 --use-conda`

Owner: Johannes Gawron
