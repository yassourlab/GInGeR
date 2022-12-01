# About GInGeR:

GInGeR is a tool for analyzing the genomic regions of genes of interest in metagenomic samples given metagenomic
sequencing data (paired-end short reads and optionally also Oxford Nanopore long reads) and a set of genes of interest.
GInGeR constructs and assembly graph, locates the genes of interest in the assembly graph, identifies their potential
genomic contexts in the graph and uses a reference genomes database to verify the contexts and assign them to carrier
species.

# Installation using conda:

0. In case you don't have conda installed, [follow the directions](https://docs.anaconda.com/anaconda/install/index.html) to install conda or miniconda (recommended).
1. Clone the repository using `git clone`.
2. Use the conda environment file given in the repository to create GInGeR's conda
   environment `conda env create -f ginger.yml`.
3. To verify your installation run the test `test_ginger.sh`.

# Running GInGeR
## supported sequencing data types

GInGer requires Illumina **paired-end** sequencing data in a fastq format (zipped fastq files are also excepted).
In addition to the short paired-end reads, it excepts long Oxford NanoPore reads (optional).

## GInGeR's command line application
To run GInGeR, you can type:
``` bash
python3 ginger.py <SHORT_READS_1> <SHORT_READS_2> <GENES_PATH>
```
where `<SHORT_READS_1>` and `<SHORT_READS_2>` are fastq (or fastq.gz) files of sequencing data and `<GENES_PATH>` is the path to the fastq of genes of interest.

In order to see additional options, `type python3 ginger.py --help` or check out our manual.

# GInGer's outputs
GInGeR's outputs two results files:
1. context_level_matches.csv - the locations in the reference database that match the genomic contexts that were
   detected in the assembly graph for the genes of interest. 
   TODO - write in detail about the format
2. species_level_matches.csv - an aggregation of the first output by gene and species, intending to summarize whether a
   gene was detected on a certain species in the sample.
   TODO - write in detail about the format

# Dependencies and requirements

# Reference database




