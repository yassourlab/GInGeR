# About GInGeR:

GInGeR is a tool for analyzing the genomic regions of genes of interest in metagenomic samples (paired-end short reads and optionally also Oxford Nanopore long reads) and a set of genes of interest (provided as **protein sequences**).
To do so, it constructs an assembly graph, locates the genes of interest in the assembly graph, identifies their potential
genomic contexts in the graph, verifies the contexts and assigns them to carrier species using a reference genomes database.

**Note that a paper describing GInGeR was not yet published. If you're interested in using it for your research, please contact netta.barak [at] mail.huji.ac.il**

# Installation using conda:

0. In case you don't have conda/miniconda/mamba
   installed, follow the directions to install [miniconda](https://docs.anaconda.com/anaconda/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).
1. Clone the repository using `git clone`.
2. Use the conda environment file given in the repository to create and activate GInGeR's conda environment:
    * `cd GInGeR`
    * `conda env create -f ginger.yml`
    * `conda activate ginger_env`
3. Install the ginger package on your conda env:
    * `python -m pip install .` (note that the `.` is part of the command)
5. Download the Kraken2 database to the GInGeR directory:
    * `wget -r -np -nH --cut-dirs=5 -R "index.html*" https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/kraken2_db_uhgg_v2.0.2/`
6. To verify your installation run ginger on a test dataset:
    * `run_ginger tests/test_files/ecoli_1K_1.fq.gz tests/test_files/ecoli_1K_2.fq.gz tests/test_files/test_gene.faa e2e_test_output --max-species-representatives 1`
    * **Note that due to Kraken2's memory requirements, you'd need to allocate at least 16G of memory for the pipeline
      to run successfully**. In case you would like to test ginger but skip the step using Kraken, you can
      run: `run_ginger tests/test_files/ecoli_1K_1.fq.gz tests/test_files/ecoli_1K_2.fq.gz tests/test_files/test_gene.faa e2e_test_output --sample-specific-references tests/test_files/merged_filtered_ref_db.fasta.gz --max-species-representatives 1`

# Running GInGeR

## Supported sequencing data types

GInGer requires Illumina **paired-end** sequencing data in fastq format (zipped fastq files are also accepted). In
addition to the short paired-end reads, it accepts long Oxford Nanopore reads (optional).

## GInGeR's command line application

To run GInGeR run:

``` bash
run_ginger <SHORT_READS_1> <SHORT_READS_2> <GENES_PATH> <OUT_DIR>
```

where `<SHORT_READS_1>` and `<SHORT_READS_2>` are fastq (or fastq.gz) files of sequencing data, `<GENES_PATH>` is the
path to a **FASTA file containing protein sequences (amino acid sequences, not nucleotide sequences)** of genes of interest, and `<OUT_DIR>` is the path specifying where to save the results and
additional files.

### Optional inputs:

`--long-reads` - A fastq or fastq.gzip file of Oxford Nanopore reads.

### Run parameters:

`-t` or `--threads` - An integer specifying the number of threads that will be used for running Kraken2, SPAdes and
Minimap2. Default is 1.

`--species-inclusion-threshold` - A float in the range [0,1] specifying the minimal fraction of reads that need to be mapped to a certain
species for it to be included in the analysis. This represents read-based abundance (proportion of DNA sequences in the sample), 
not cell count abundance. Species with larger genomes contribute more reads than those with smaller genomes at equal cell counts. 
Default is 0.01 (1%).

`--max-species-representatives` - The maximal references per species that will be downloaded from UHGG and taken into
account in the aggregation of results at the species level. Default is 100.

`--depth-limit` - An integer specifying the maximal depth for paths describing context candidates in the assembly graph.
Default is 12. Greater values can increase runtime and generate longer and lower certainty contexts. Smaller values can
decrease runtime and generate shorter contexts with higher certainty.

`--max-context-len` - The required length of one-sided context candidate. Default is 2500. Shorter contexts can be
generated if the algorithm reaches a dead end in the assembly graph or if the path reaches the `--depth-limit`.

`--gene-pident-filtering-th` - A float in the range [0,1] specifying the minimal % of matched base pairs required for
locating a gene in the graph. Default is 0.9.

`--paths-pident-filtering-th` - A float in the range [0,1] specifying the minimal % of matched base pairs required for
matching a context candidate to a reference sequence. Default is 0.9.

`--return-all-gene-matches` - A flag indicating whether to return all gene matches without applying non-max-suppression (NMS). 
By default (when this flag is not set), GInGeR applies NMS to filter overlapping gene matches, keeping only the top-scoring 
matches. When this flag is set, all gene matches above the `--gene-pident-filtering-th` threshold are returned.

`--nms-iou-threshold` - A float in the range [0,1] specifying the IoU (Intersection over Union) threshold for 
non-max-suppression when filtering overlapping gene matches. Gene matches with IoU > this threshold are considered 
overlapping, and only the highest-scoring match is kept. Only used when `--return-all-gene-matches` is not set. 
Default is 0.8.

### Advanced options (skipping steps or customizing reference database):

`--skip-assembly` - A flag that indicates whether to skip the assembly step. If the flag is set, the
argument `--assembly-dir` must be supplied and point to the results of a SPAdes run.

`--assembly-dir` - Specifies where to save the assembly results. Default is `<OUT_DIR>/SPAdes`. In case you want to use
a pre-ran assembly, please specify here the directory of SPAdes' output.

`--kraken-output-path` - A path for saving Kraken2's output..

`--reference-genomes-metadata` - The path to the reference database metadata table (see more information in the Reference database
section).

`--downloaded-references-dir` - The directory to which GInGeR will download missing reference genomes from UHGG. Defaults
to `references_dir`. This folder can be shared for all runs of GInGeR in order to avoid the same file being downloaded
and saved multiple times.

`--sample-specific-references` - A reference database in a fasta format. Using this will skip the stages of creating a sample-specific database based on the species Kraken2 detected in the sample.

`--keep-intermediate` - Controls which intermediate files are kept after the pipeline completes. By default, all files are kept (`all`). Options include: `final` (only result CSVs), `assembly` (SPAdes output), `alignment` (PAF/M8 files), `sequences` (FASTA files), `kraken` (Kraken2/Bracken output), `reference` (reference database files), `plasmid` (GeNomad input fasta and output directory). Can specify multiple categories by repeating the flag. Example: `--keep-intermediate final --keep-intermediate assembly` keeps only the final CSV results and the SPAdes directory.

### Plasmid scoring

`--add-plasmid-score` / `--no-add-plasmid-score` - Whether to run [GeNomad](https://github.com/apcamargo/genomad)
on the genomic contexts found for each gene (and on the contigs of genes that were located in the
assembly but had no species-level match) and add `plasmid_score` columns to the output CSVs. Default: enabled.

`--genomad-db` - The path to GeNomad's database directory. Required when `--add-plasmid-score` is set (the
default). Create it with `genomad download-database <path>` (see
[GeNomad's documentation](https://portal.nersc.gov/genomad/)). Defaults to a `genomad_db` directory next to the
GInGeR package.

# GInGeR's outputs

GInGeR outputs the following result files:

1. **context_level_matches.csv** - A CSV file specifying the locations in the reference database that match the genomic contexts
   that were detected in the assembly graph for the genes of interest. The CSV columns are:
    * gene - the gene of interest (name is taken from the fasta of genes given as input to the pipeline)
    * reference_contig - the reference sequence (usually an assembly contig) where the context of gene in the sample was
      matched
    * in_context - the name of the incoming context that was matched (the equivalent sequence can be found in
      all_in_paths.fasta)
    * out_context - the name of the outgoing context that was matched (the equivalent sequence can be found in
      all_out_paths.fasta)
    * in_context_start - the start of the match of the incoming context to the reference sequence
    * gene_start - the end of the match of the incoming context to the reference sequence (note that the gene is not
      necessarily found there, this is the estimated start location of the gene)
    * gene_end - the start of the match of the outgoing context to the reference sequence
    * out_context_end - the end of the match of the outgoing context to the reference sequence
    * match_score - the match score. Namely, the % of matching base-pairs between the genomic contexts and the reference
      sequence
    * genome - the genome id as inferred from the reference_contig field
    * species - the species name as inferred from UHGG-metadata.tsv
    * plasmid_score - (only when `--add-plasmid-score` is set, default) [GeNomad](https://github.com/apcamargo/genomad)'s
      plasmid score (in the range [0,1]) for the sequence formed by the in_context, the gene and the out_context. 0 if
      GeNomad did not report a score for this sequence.
    * context_species_diversity - the Shannon diversity index of the species matched by this row's in_context/out_context
      pair (the in-gene-out trio). For each species matched by the trio, the number of unique reference genomes it was
      matched to is divided by that species' number of references (capped at `--max-species-representatives`), and the
      resulting corrected counts are normalized to probabilities before computing the Shannon diversity. A trio matched
      to references of a single species has a diversity of 0; the more evenly a trio is spread across distinct species,
      the higher its diversity.
    * context_species_confidence_score - how confidently this row's in_context/out_context pair (the in-gene-out trio)
      points to this row's species, versus other species also matched by the same trio. in_context and out_context are
      scored independently using the same logic: for every species matched by the context, the number of matches is
      divided by that species' number of references (capped at `--max-species-representatives`), and the resulting
      corrected counts are normalized so the scores of all species matched by that context sum to 1.
      context_species_confidence_score is the average of the in_context and out_context scores.

2. **species_level_matches.csv** - An aggregation of the first output by gene and species, summarizing whether a
   gene was detected in a certain species in the sample. The CSV columns are:
    * gene - the gene of interest as described for `context_level_output.csv`
    * species - the species name as described for `context_level_output.csv`
    * references_ratio - the % of instances of the species that included the gene (given as a ratio in the range [0,1])
    * match_score_max - the maximal match score for the given gene and species
    * species_instances - the number of instances of the given species taken into account in the calculation (determined
      as the minimum between `--max-species-representatives` and the number of instances found in the UHGG database)
    * plasmid_score_mean - (only when `--add-plasmid-score` is set, default) the average `plasmid_score` across all of
      the gene's distinct genomic contexts matched to this species
    * plasmid_score_most_common_context - (only when `--add-plasmid-score` is set, default) the `plasmid_score` of the
      genomic context matched to the largest number of reference genomes of this species (averaged across ties)
    * n_contexts - the number of distinct in-gene-out trios (in_context/out_context pairs) that contributed to this
      gene-species match
    * species_confidence_score - the mean of `context_species_confidence_score` over this gene-species pair's distinct
      in-gene-out trios (in_context/out_context pairs)

3. **genes_detected_in_graph_with_no_species_match.csv** - Genes that were located in the assembly graph but for which
   no genomic context could be matched to a reference genome (and so don't appear in the outputs above). The CSV
   columns are:
    * gene - the gene of interest as described for `context_level_output.csv`
    * contig - the assembly contig the gene was found on
    * gene_match_score - the gene-to-contig match score (the % of matching base-pairs between the gene and the contig)
    * plasmid_score - (only when `--add-plasmid-score` is set, default) GeNomad's plasmid score for the full contig the
      gene was found on. 0 if GeNomad did not report a score for this contig.

# Reference database

By default, GInGeR uses the Unified Human Gastrointestinal Genome collection (UHGG,
see [genome catalog](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/README_v2.0.txt) and
[Almeida et al.](https://www.nature.com/articles/s41587-020-0603-3)) as a reference database. When
running GInGer, it first detects the dominant species in the sample using Kraken2 and estimates per-species read counts using Bracken. It then keeps only species with estimated sequencing coverage above `--species-coverage-threshold` (default: 10). Estimated coverage is computed as:

new_est_reads * (avg_len_R1 + avg_len_R2) / median_genome_length

where `new_est_reads` comes from Bracken, read lengths come from `seqkit stats`, and `median_genome_length` is computed from the top UHGG reference genomes per species (ranked by the metadata Quality score) capped by `--max-species-representatives`. GInGeR then downloads missing
references (`--max-species-representatives` of references per species) for these species from UHGG to the `--downloaded-references-dir` and
combining them into a reference database customized for the given sample. It is recommended to use a
shared `--downloaded-references-dir` for multiple samples (this is GInGeR's default behavior) to save time and storage by
not downloading the same references multiple times.

### Supply GInGeR with a FASTA of reference species

You can do so using the `--sample-specific-references` option. For optimal results, it should include only species that
are relevant for your sample. In this case you should also supply a `--reference-genomes-metadata` file in TSV format with the header
'Genome Length Lineage FTP_download species':

- 'Genomes' should be the id of the genome as appears in your fasta file. Each contig should have the following format
  <genome_id>_<contig_num>.
- 'Completeness' can be left with Null values because it is not used in this option.
- 'FTP_download' can be left with Null values because it is not used in this option.
- 'species' the species name of the given genome. Will be used when generating the outputs.

# Dependencies

GInGeR requires a 64-bit Linux system and conda (or miniconda).

### Python 3.12.13

The pipeline might also work on other versions of Python, but they were not tested.

#### Python packages

  - biopython=1.87
  - pandas=3.0.3
  - kraken2=2.17.1
  - bracken=3.1
  - minimap2=2.30
  - pip=26.1.1
  - spades=4.2.0
  - click=8.4.0
  - pyfastg=0.2.0
  - pafpy=0.2.0
  - tqdm=4.67.3
  - mmseqs2=18.8cc5c
  - genomad=1.12.0
  - seqkit=2.13.0

# Acknowledgements

GInGeR is open source and available thanks to the hard work of [Haimasree Bhattacharya](https://github.com/haimasree).




