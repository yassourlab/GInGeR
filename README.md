# About GInGeR:

GInGeR is a tool for analyzing the genomic regions of genes of interest in metagenomic samples given metagenomic
sequencing data (paired-end short reads and optionally also Oxford Nanopore long reads) and a set of genes of interest.
GInGeR constructs and assembly graph, locates the genes of interest in the assembly graph, identifies their potential
genomic contexts in the graph and uses a reference genomes database to verify the contexts and assign them to carrier
species.

# Installation using conda:

0. In case you don't have conda
   installed, [follow the directions](https://docs.anaconda.com/anaconda/install/index.html) to install conda or
   miniconda (recommended).
1. Clone the repository using `git clone`.
2. Use the conda environment file given in the repository to create GInGeR's conda
   environment `conda env create -f ginger.yml`.
3. To verify your installation run the test `test_ginger.sh`.

# Running GInGeR

## supported sequencing data types

GInGer requires Illumina **paired-end** sequencing data in a fastq format (zipped fastq files are also excepted). In
addition to the short paired-end reads, it excepts long Oxford NanoPore reads (optional).

## GInGeR's command line application

To run GInGeR, you can type:

``` bash
python3 ginger.py <SHORT_READS_1> <SHORT_READS_2> <GENES_PATH> <OUT_DIR>
```

where `<SHORT_READS_1>` and `<SHORT_READS_2>` are fastq (or fastq.gz) files of sequencing data, `<GENES_PATH>` is the
path to the fastq of genes of interest, and `<OUT_DIR>` is the path specifying the where to save the results and
additional files.

### Optional inputs:

`--long-reads` - A fastq or fastq.gzip file of Oxford Nanopore reads.

### Run parameters:

`-t` or `--threads` - An integer specifying the Number of threads that will be used for running Kraken2, SPAdes and
Minimap2.

`--reads-ratio-th` - A float in the range [0,1] specifying the minimal % of reads that need to be mapped to a certain
species for it to be included in the analysis, default 0.01 (1%).

`--depth-limit` - An integer specifying The maximal depth for paths describing context candidates in the assembly graph.
default is 12. Greater values can increase runtime and generate longer and lower certainty contexts. Smaller values can
decrease runtime and generate shorter contexts with higher certainty.

`--max-context-len` - The required length of one-sided context candidate. Default is 2500. Shorter contexts can be
generated if the algorithm reaches a dead end in the assembly graph or if the path reaches the `--depth-limit`.

`--gene-pident-filtering-th` - A float in the range [0,1] specifying the minimal % of matched base pairs required for
locating a gene in the graph

`--paths-pident-filtering-th` - A float in the range [0,1] specifying the minimal % of matched base pairs required for
matching a context candidate to a reference sequence

### Advanced options (skipping steps or customizing reference database):

`--skip-assembly` - A flag that indicates whether to skip the assembly step. If the flag is set to True, the
argument `--assembly--dir` must be supplied and direct to the results of a SPAdes run,

`--assembly-dir` - Specifies where to save the assembly results, default is `<OUT_DIR>/SPAdes`. In case you want to use
a pre-ran assembly, please specify here the directory of sapdes' output.

`--kraken-output-path` - A path for saving Kraken2's output

`--metadata-path` - The path to the reference database metadata table. (see more information in the Reference database
section)

`--references-dir`- The directory to which GInGeR will download missing reference genomes from UHGG. Defaults
to `references_dir`. This folder can be shared for all runs of GInGer in order to avoid the same file being downloaded
and saved multiple times.

`--merged-filtered-fasta` - A reference database in a fasta format. using this will skip the stages of creating a sample
specific database based on the species Kraken2 detected in the sample.

# GInGer's outputs

GInGeR's outputs two results files:

1. context_level_matches.csv - the locations in the reference database that match the genomic contexts that were
   detected in the assembly graph for the genes of interest. TODO - write in detail about the format
2. species_level_matches.csv - an aggregation of the first output by gene and species, intending to summarize whether a
   gene was detected on a certain species in the sample. TODO - write in detail about the format

# Reference database

By default, GInGeR uses the Unified Human Gastrointestinal Genome collection (UHGG,
see [genome catalog](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/README_v2.0.txt) and
[Publication by Almeida et al.](https://www.nature.com/articles/s41587-020-0603-3)) as a reference database. When
running GInGer, it first detects the dominant species in the sample using Kraken2, and then downloads missing
references (a limited number of references per species) for these species from UHGG to the `--references-dir` and
combining them into a reference database customized for the given sample. It is recommended to use a
shared `--references-dir` for multiple samples (this is GInGeR's default behavior) in order to save time and storage by
not downloading the same references multiple times.

## Using a different reference database

UHGG is a comprehensive and high quality database that includes over 200K assemblies of ~4.7K species found in the human
gut. Therefore, if you work with Human gut microbiome samples, UHGG is probably a great choice of a reference database.
In case you're working with samples from other environments, you might want to run GInGeR using a different reference
database. In order to do so, you can choose one of two options (in both cases, you would need to adapt your database to
match UHGG's style):

### supply GInGer with a fasta of reference species

You can do so using using the `--merged-filtered-fasta` option. For optimal results, it should include only species that
are relevant for your sample. In this case you should supply also a `--metadata-path` in a TSV format and the header '
Genome Length Lineage FTP_download species':

- 'Genomes' should be the id of the genome as appears in your fasta file. Each contig should have the following format
  <genome_id>_<contig_num>.
- 'Length' can be left with Null values because it is not used in this option.
- 'Lineage' can be left with Null values because it is not used in this option.
- 'FTP_download' can be left with Null values because it is not used in this option.
- 'species' the species name of the given genome. Will be used when generating the outputs.

### Use GInGer's pipeline to filter and select the relevant species from your reference database
TODO - missing this part
# Dependencies (should all be automatically included in the conda environment)

### Python 3.10

The pipeline might also work on older versions of python, but they were not tested.

#### python packages

- biopython=1.78
- pandas=1.4.4
- click=8.0.4
- pyfastg==0.1.0
- readpaf==0.0.10

### SPAdes - version 3.15.5

Used for building an assembly graph for the given metagenomic sample.

### Kraken2 - version 2.1.2

Used for selecting only the references of the species that are found in the sample.

### Minimap2 - version 2.24

Used for searching for locating the genes in the assembly graphs (by searching for them in the assembly contigs and
inferring their location in the graph)






