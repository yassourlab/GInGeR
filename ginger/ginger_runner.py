import sys
import logging
import click
import os
import pandas as pd
import shutil
from glob import glob

from ginger import locating_genes_in_graph as lg
from ginger import reference_database_utils as rdu
from ginger import assembly_utils as au
from ginger import extract_contexts_candidates as ecc
from ginger import sequence_alignment_utils as sau
from ginger import verify_context_candidates as vcc
from ginger import pipeline_utils as pu
from ginger import plasmid_detection_utils as pdu
from ginger import constants as c


logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

def resolve_reference_source(reference_source, kraken_db, reference_genomes_metadata,
                             downloaded_references_dir):
    """Resolve the settings that go with a reference catalog.

    The three paths are filled in only where the caller left them unset - the metadata table ships
    inside the package and the Kraken database sits beside it in the repo, which is where the
    previous per-option defaults pointed, while the references directory stays relative to the
    working directory. The distinct-k-mer ratio threshold always comes from the catalog, since the
    right value depends on how finely that catalog splits species.
    """
    source = rdu.REFERENCE_SOURCES[reference_source]
    package_dir = os.path.dirname(__file__)
    if kraken_db is None:
        kraken_db = os.path.join(package_dir, '..', source['kraken_db'])
    if reference_genomes_metadata is None:
        reference_genomes_metadata = os.path.join(package_dir, source['metadata'])
    if downloaded_references_dir is None:
        downloaded_references_dir = source['references_dir']
    distinct_kmer_ratio_threshold = source['distinct_kmer_ratio_threshold']
    log.info(f'reference source {reference_source}: metadata {reference_genomes_metadata}, '
             f'kraken db {kraken_db}, references dir {downloaded_references_dir}, '
             f'distinct k-mer ratio threshold {distinct_kmer_ratio_threshold}')
    return (kraken_db, reference_genomes_metadata, downloaded_references_dir,
            distinct_kmer_ratio_threshold)


def cleanup_intermediate_files(out_dir, keep_options):
    """Remove intermediate files based on keep_options."""
    if 'all' in keep_options:
        return
    
    log.info('Cleaning up intermediate files')
    
    # Define file patterns for each category.
    # in_gene_out_contexts.fasta is deliberately in none of them - it is a result, not an
    # intermediate, and context_level_matches.csv is unusable for visualization without it.
    # --write-context-sequences is what decides whether it is there at all
    cleanup_map = {
        'assembly': ['SPAdes'],
        'alignment': ['*.paf', '*.m8', 'mmseqs_tmp', 'nodes_to_contigs_w_gaps.paf'],
        'sequences': ['all_in_paths.fasta', 'all_out_paths.fasta'],
        'kraken': ['kraken_*.tsv', 'bracken_*.tsv'],
        'reference': ['merged_filtered_ref_db.*', 'references_used.csv', 'reference_contig_to_genome.tsv'],
        # genomad_output only survives a GeNomad that failed - a successful run keeps just the summary
        'plasmid': ['plasmid_summary.tsv', 'genomad_output'],
    }
    
    # Remove categories not in keep_options
    # 'final' is not a real category, it just means "keep final CSV results"
    for category, patterns in cleanup_map.items():
        if category in keep_options:
            continue
        for pattern in patterns:
            # glob handles an exact name too - it yields the path when it exists and nothing when it does not
            for path in glob(os.path.join(out_dir, pattern)):
                try:
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                        log.debug(f'Removed directory {path}')
                    else:
                        os.remove(path)
                        log.debug(f'Removed file {path}')
                except Exception as e:
                    log.warning(f'Failed to remove {path}: {e}')


@click.command()
@click.argument('short-reads-1', required=True, type=click.Path(exists=True))
@click.argument('short-reads-2', required=True, type=click.Path(exists=True))
@click.argument('genes-path', required=True, type=click.Path(exists=True))
@click.argument('out-dir', required=True, type=click.Path())
@click.option('--long-reads', type=click.Path(exists=True), help='A fastq or fastq.gzip file of Oxford Nanopore reads')
@click.option('--assembly-dir', default=None,
              help="Specifies where to save the assembly results. In case of pre-ran assembly, please insert the path do the spades output directory")
@click.option('--threads', '-t', type=int, default=1,
              help='Number of threads that will be used for running Kraken2, SPAdes and Minimap2')
@click.option('--kraken-report-path', type=click.Path(exists=True), default=None,
              help="A Kraken2 report from a previous run on the same reads, which skips running Kraken2 (Bracken and the later stages still run). It must have been created with --report-minimizer-data, and with the same database as --kraken-db")
@click.option('--reference-source', type=click.Choice(sorted(rdu.REFERENCE_SOURCES)),
              default=rdu.DEFAULT_REFERENCE_SOURCE,
              help="Which reference catalog to use. Selects how reference genomes are downloaded - GTDB's by assembly accession from NCBI, UHGG's as .gff.gz over FTP - and supplies matching defaults for --kraken-db, --reference-genomes-metadata and --downloaded-references-dir. Any of those passed explicitly wins.")
@click.option('--kraken-db', type=click.Path(), default=None,
              help='The path to the Kraken2 database directory. Defaults to the one matching --reference-source')
@click.option('--species-coverage-threshold', type=float, default=3,
              help='The minimal estimated sequencing coverage required for including a species in the analysis. Coverage is estimated as: bracken_estimated_reads * (avg_len_R1 + avg_len_R2) / median_genome_length, where median genome length is computed from the top references per species (by Quality) capped by --max-species-representatives. Default 3.')
@click.option('--max-species-representatives', type=int, default=100,
              help='The maximal references per species that will be downloaded from NCBI and taken into account in the aggregation of results at the species level')
@click.option('--reference-genomes-metadata', type=click.Path(), default=None,
              help='The path to the reference database metadata table. Defaults to the one matching --reference-source')
@click.option('--downloaded-references-dir', type=click.Path(), default=None,
              help="The directory to which GInGeR will download missing reference genomes. This folder can be shared for all runs of GInGer in order to avoid the same file being  downloaded and saved multiple times. Defaults to references_dir_{--reference-source}, so the two catalogs never share a cache")
@click.option('--sample-specific-references', type=click.Path(), default=None,
              help='A fasta, fasta.gz or mmi (minimap indexed) file that will be used a reference database (using this will skip the stages of creating a sample specific database based on the species detected in the sample by Kraken2)')
@click.option('--reference-contig-to-genome', type=click.Path(), default=None,
              help="A TSV with 'contig' and 'Genome' columns mapping every contig of --sample-specific-references to the reference genome it belongs to. GInGeR writes one (reference_contig_to_genome.tsv) whenever it builds the reference database itself, so pass that file back when reusing a database across runs. Without it the genome is guessed from the contig name, which only works when contigs are named {genome_id}_{contig_num}")
@click.option('--depth-limit', type=int, default=12,
              help='The maximal depth for paths describing context candidates in the assembly graph')
@click.option('--max-gap-ratio', type=float, default=1.5,
              help="The maximal ratio between the length of the gene and the gap between it's contexts in the database")
@click.option('--context-len', type=int, default=2500,
              help='The length of a one-sided context candidate. Contexts are always exactly this long - a side that cannot supply that much sequence gets no context, and a gene needs a context on both sides to be reported.')
@click.option('--gene-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for locating a gene in the graph')
@click.option('--paths-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for matching a context candidate to a reference sequence')
@click.option('--keep-intermediate', multiple=True,
              type=click.Choice(['all', 'final', 'assembly', 'alignment', 'sequences', 'kraken', 'reference', 'plasmid'],
                               case_sensitive=False),
              default=['all'],
              help="Specify which intermediate files to keep. Options: all (default, keep everything), final (only result files), assembly (SPAdes output), alignment (PAF/M8 files), sequences (FASTA files), kraken (Kraken2/Bracken output), reference (reference database files), plasmid (GeNomad's plasmid summary). Can specify multiple by repeating the flag: --keep-intermediate final --keep-intermediate assembly. in_gene_out_contexts.fasta is a result, so when it is written (see --write-context-sequences) it is always kept.")
@click.option('--skip-assembly', is_flag=True, default=False,
              help='A flag that indicates whether or not to skip the assembly step. If the flag is set to True, the argument --assembly--dir must be supplied and direct to the results of a SPAdes run')
@click.option('--return-all-gene-matches', is_flag=True, default=False,
              help='By default, GInGeR applies non-max-suppression (NMS) to the alignment of genes of interest to the assembly graph, keeping only top scoring matches and removing redundant overlapping matches. If this flag is set to True, all gene matches will be returned without applying NMS.')
@click.option('--nms-iou-threshold', type=float, default=0.8,
              help='The IoU (Intersection over Union) threshold used for non-max-suppression when filtering overlapping gene matches. Gene matches with IoU > this threshold are considered overlapping. Only used when --return-all-gene-matches is False. Default: 0.8')
@click.option('--add-plasmid-score/--no-add-plasmid-score', default=True,
              help='Run GeNomad on the genomic contexts found for each gene (and on the contigs of genes with no species-level match) and add plasmid_score columns to the output CSVs. Requires --genomad-db to point to a valid GeNomad database. Default: True')
@click.option('--genomad-db', type=click.Path(),
              default=os.path.join(os.path.dirname(__file__), '..', 'genomad_db'),
              help="The path to GeNomad's database directory (create one with `genomad download-database <path>`). Only used when --add-plasmid-score is set.")
@click.option('--contig-context-fallback/--no-contig-context-fallback', default=True,
              help='For a gene found on a gap-containing contig (a contig SPAdes assembled from several graph paths joined using paired-end evidence), also take its context from the flanking sequence of the contig itself. The assembly graph describes such a gene\'s context poorly or not at all, but a context taken from the contig may cross one of those joins rather than a graph edge, so it is named "{gene}|{contig}|{start}|{end}|{score}|{nodes}|contigfallback|{side}" in the output. Default: True')
@click.option('--write-context-sequences/--no-write-context-sequences', default=False,
              help='Write in_gene_out_contexts.fasta - the in-gene-out sequence behind every context level row, and the full contig of every gene with no species-level match - along with the context_seq_id and gene offset columns that join a row to its sequence. Tens of MB for a typical sample, so it is off unless asked for. Implied by --add-plasmid-score, which needs this fasta as GeNomad\'s input. Default: False')
def run_ginger_e2e(**kwargs):
    """GInGeR - A tool for analyzing the genomic contexts of genes in metagenomic samples.

    \b
    SHORT_READS_1 - R1 fastq or fastq.gzip file

    \b
    SHORT_READS_2 - R2 fastq or fastq.gzip file

    \b
    GENES_PATH - A fasta file with the genes of interest

    \b
    OUT_DIR - A path specifying where to save GInGeR's output

    """
    # click passes every argument and option by keyword, under the same names ginger_e2e_func takes
    return ginger_e2e_func(**kwargs)


def ginger_e2e_func(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads, kraken_report_path,
                    kraken_db, species_coverage_threshold, reference_genomes_metadata, downloaded_references_dir, sample_specific_references, genes_path, depth_limit,
                    max_gap_ratio, context_len, gene_pident_filtering_th,
                    paths_pident_filtering_th, keep_intermediate, skip_assembly, max_species_representatives, return_all_gene_matches, nms_iou_threshold,
                    add_plasmid_score=True, genomad_db=None, contig_context_fallback=True,
                    write_context_sequences=False, reference_contig_to_genome=None,
                    reference_source=rdu.DEFAULT_REFERENCE_SOURCE):
    # Log the command that was run
    log.info(f"Running GInGeR with command: {' '.join(sys.argv)}")

    # whatever the caller did not pin down comes from the chosen catalog, so a GTDB metadata table
    # can't end up paired with the UHGG Kraken database by omission
    kraken_db, reference_genomes_metadata, downloaded_references_dir, distinct_kmer_ratio_threshold = \
        resolve_reference_source(reference_source, kraken_db, reference_genomes_metadata,
                                 downloaded_references_dir)

    pu.ensure_out_dir_is_fresh(out_dir)
    # create output directory if it doesn't exist
    pu.check_and_make_dir_no_file_name(out_dir)
    # filter reference database using kraken
    references_used_path = c.REFERENCES_USED_TEMPLATE.format(out_dir=out_dir)
    # how a match on a reference contig is traced back to a genome and so to a species. GInGeR
    # writes it when it builds the reference database; with --sample-specific-references the user
    # supplies the one saved from the run that built that database, or leaves it unset and lets the
    # genome be guessed from the contig name
    contig_to_genome_path = reference_contig_to_genome
    if sample_specific_references is None:
        sample_specific_references = c.MERGED_FILTERED_REF_DB_TEMPLATE.format(out_dir=out_dir)
        contig_to_genome_path = c.CONTIG_TO_GENOME_TEMPLATE.format(out_dir=out_dir)
        kraken_output_path = c.KRAKEN_OUTPUT_TEMPLATE.format(out_dir=out_dir)
        existing_kraken_report = kraken_report_path
        kraken_report_path = c.KRAKEN_REPORT_TEMPLATE.format(out_dir=out_dir)
        bracken_output = c.BRACKEN_OUTPUT_TEMPLATE.format(out_dir=out_dir)
        bracken_report = c.BRACKEN_REPORT_TEMPLATE.format(out_dir=out_dir)
        species_included_in_analysis_path = c.SPECIES_INCLUDED_IN_ANALYSIS_TEMPLATE.format(out_dir=out_dir)
        rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path,
                                             kraken_report_path, bracken_output, bracken_report, species_coverage_threshold,
                                             reference_genomes_metadata, downloaded_references_dir, sample_specific_references,
                                             references_used_path,
                                             max_species_representatives, kraken_db,
                                             species_included_in_analysis_path, contig_to_genome_path,
                                             existing_kraken_report=existing_kraken_report,
                                             source=reference_source,
                                             distinct_kmer_ratio_threshold=distinct_kmer_ratio_threshold)
    if not sample_specific_references.endswith('mmi'):
        indexed_reference = sau.generate_index(sample_specific_references, sau.INDEXING_PRESET)
    else:
        indexed_reference = sample_specific_references
    # run assembly
    if assembly_dir is None:
        assembly_dir = c.ASSEMBLY_DIR_TEMPLATE.format(out_dir=out_dir)
    if not skip_assembly:
        au.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_dir, threads)
    # run tool

    assembly_graph, genes_to_analyze, geometry, contigs_with_gaps = lg.locate_genes_in_graph(assembly_dir,
                                                                                                  gene_pident_filtering_th,
                                                                                                  genes_path,
                                                                                                  threads,
                                                                                                  out_dir,
                                                                                                  return_all_gene_matches,
                                                                                                  nms_iou_threshold)
    if not genes_to_analyze:
        log.info(
            'No genes of interest detected in assembly. GInGeR run stopped - no results generated')
        return

    genes_detected_no_species_match_output_path = c.GENES_DETECTED_IN_GRAPH_WITH_NO_SPECIES_MATCH_OUTPUT_TEMPLATE.format(
        out_dir=out_dir
    )

    # get in and out paths
    in_paths_fasta = c.IN_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    out_paths_fasta = c.OUT_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(
                                                              assembly_graph, geometry,
                                                              genes_to_analyze, depth_limit,
                                                              context_len, in_paths_fasta,
                                                              out_paths_fasta,
                                                              c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=assembly_dir),
                                                              contigs_with_gaps if contig_context_fallback else frozenset())

    # map them to the reference
    in_contexts_to_ref_genomes = c.IN_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE.format(temp_folder=out_dir)
    out_contexts_to_ref_genomes = c.OUT_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE.format(temp_folder=out_dir)

    sau.map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, indexed_reference, in_contexts_to_ref_genomes,
                                       out_contexts_to_ref_genomes, threads)

    # merge and get results
    # built once and handed to both stages, so that the species a match was resolved to and the
    # Genome column written beside it in the CSV always come from the same mapping
    contig_species_lookup = vcc.build_contig_species_lookup(reference_genomes_metadata, contig_to_genome_path)
    context_level_results = vcc.process_in_and_out_paths_to_results(in_contexts_to_ref_genomes,
                                                                    out_contexts_to_ref_genomes,
                                                                    gene_lengths, paths_pident_filtering_th, 0,
                                                                    max_gap_ratio, reference_genomes_metadata,
                                                                    contig_species_lookup=contig_species_lookup)

    # write the sequence behind every conclusion below - the in-gene-out sequence of every context
    # level row, and the contig of every gene with no context match. Tens of MB for a typical sample,
    # so it is written when the user asked for it to visualize and analyze, and when GeNomad is going
    # to run on it either way
    contexts_fasta_path, context_seq_records = None, {}
    if write_context_sequences or add_plasmid_score:
        genes_with_context_matches = {gene for gene, _ in context_level_results.keys()} if context_level_results else set()
        contexts_fasta_path, context_seq_records = pu.write_in_gene_out_contexts_fasta(
            context_level_results, genes_to_analyze, genes_with_context_matches,
            in_paths_fasta, out_paths_fasta, c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=assembly_dir),
            c.IN_GENE_OUT_CONTEXTS_FASTA_TEMPLATE.format(out_dir=out_dir))

    # run GeNomad on the gene contexts and on the contigs of genes with no species-level match
    context_plasmid_scores, contig_plasmid_scores = None, None
    if add_plasmid_score and contexts_fasta_path:
        genomad_out_dir = c.GENOMAD_OUTPUT_DIR_TEMPLATE.format(out_dir=out_dir)
        plasmid_summary_path = pdu.run_genomad(contexts_fasta_path, genomad_out_dir, genomad_db, threads)
        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(plasmid_summary_path,
                                                                                 context_seq_records)
        pdu.keep_only_plasmid_summary(genomad_out_dir, plasmid_summary_path,
                                      c.PLASMID_SUMMARY_TEMPLATE.format(out_dir=out_dir))

    if not context_level_results:
        wrote_no_species_match_csv = pu.write_genes_detected_in_graph_with_no_species_match(
            genes_to_analyze,
            matched_genes=set(),
            csv_path=genes_detected_no_species_match_output_path,
        )
        if wrote_no_species_match_csv and contig_plasmid_scores is not None:
            pdu.add_plasmid_scores_to_genes_no_species_match_csv(genes_detected_no_species_match_output_path,
                                                                 contig_plasmid_scores)
        log.info(
            'No matching pairs of incoming and outgoing contexts found in reference sequences. GInGeR run stopped - no results generated')
        return
    context_level_output_path = c.CONTEXT_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    species_level_output_path = c.SPECIES_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    subspecies_level_output_path = c.SUBSPECIES_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    pu.write_context_level_output_to_csv(context_level_results, context_level_output_path, reference_genomes_metadata,
                                          max_species_representatives, contig_species_lookup=contig_species_lookup)
    # the join key onto in_gene_out_contexts.fasta, added before the plasmid scores so that a row
    # carries its sequence id even when GeNomad did not run. Only when that fasta was written - a
    # column pointing into a file that does not exist would be worse than no column
    if context_seq_records:
        pu.add_context_seq_ids_to_context_level_csv(context_level_output_path, context_seq_records)
    if context_plasmid_scores is not None:
        pdu.add_plasmid_scores_to_context_level_csv(context_level_output_path, context_plasmid_scores)

    species_level_df = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path,
                                                                                               reference_genomes_metadata,
                                                                                               species_level_output_path,
                                                                                               max_species_representatives)

    if os.path.exists(references_used_path) and 'subspecies' in pd.read_table(references_used_path).columns:
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path,
                                                                                references_used_path,
                                                                                subspecies_level_output_path,
                                                                                max_species_representatives,
                                                                                'subspecies')

    matched_genes = set()
    if species_level_df is not None and len(species_level_df.index) > 0:
        # species_level_df has a MultiIndex (gene, species)
        matched_genes = set(species_level_df.index.get_level_values(0))

    wrote_no_species_match_csv = pu.write_genes_detected_in_graph_with_no_species_match(
        genes_to_analyze,
        matched_genes=matched_genes,
        csv_path=genes_detected_no_species_match_output_path,
    )
    if wrote_no_species_match_csv and contig_plasmid_scores is not None:
        pdu.add_plasmid_scores_to_genes_no_species_match_csv(genes_detected_no_species_match_output_path,
                                                             contig_plasmid_scores)
    # Clean up intermediate files if requested
    cleanup_intermediate_files(out_dir, keep_intermediate)
    
    log.info(
        f"GInGeR completed successfully. Context-level output: {context_level_output_path}, Species-level output: {species_level_output_path}")

if __name__ == "__main__":
    run_ginger_e2e()
