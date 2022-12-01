import sys
import logging
import click
from Bio import SeqIO

import locating_genes_in_graph as lg
import reference_database_utils as rdu
import hybrid_assembly_utils as hau
import extract_contexts_candidates as ecc
import sequence_alignment_utils as sau
import verify_context_candidates as vcc
import pipeline_utils as pu
import constants as c

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


def extract_genes_lengths(genes_path):
    gene_lengths = {}
    with open(genes_path) as f:
        fasta_sequences = SeqIO.parse(f, 'fasta')
        for s in fasta_sequences:
            gene_length = len(s.seq)
            gene_name = s.id

            gene_lengths[gene_name] = gene_length
    return gene_lengths


@click.command()
@click.argument('short-reads-1', required=True, type=click.Path(exists=True))
@click.argument('short-reads-2', required=True, type=click.Path(exists=True))
@click.argument('genes-path', required=True, type=click.Path(exists=True))
@click.option('--long-reads', type=click.Path(exists=True), help='A fastq or fastq.gzip file of Oxford Nanopore reads')
@click.option('-o', '--out-dir', default='ginger_output', help="A path specifying where to save GInGeR's output",
              type=click.Path())
@click.option('--assembly-dir', default='ginger_output/SPAdes',
              help="Specifies where to save the assembly results. In case of pre-ran assembly, please insert the path do the spades output directory",
              type=click.Path())
@click.option('--threads', '-t', type=int, default=1,
              help='Number of threads that will be used for running Kraken2, SPAdes and Minimap2')
@click.option('--kraken-output-path', default=None, help="A path for saving Kraken2's output")
@click.option('--reads-ratio-th', type=float, default=0.01,
              help='The minimal % of reads that need to be mapped to a certain species for it to be included in the analysis')
@click.option('--metadata-path', type=click.Path(),
              default='UHGG-metadata.tsv')  # TODO figure out how do I give it the correct default path when I don't run GInGeR from the GInGeR dir
@click.option('--references-dir', type=click.Path(), default='ginger_output/references_folder',
              help='The directory to which GInGeR will download missing reference genomes from UHGG. This folder can be shared for all runs of GInGer in order to avoid the same file being  downloaded and saved multiple times')
@click.option('--merged-filtered-fasta', type=click.Path(), default=None,
              help='A reference database (using this will skip the stages of creating a sample specific database based on the species Kraken2) detected in the sample')
@click.option('--depth-limit', type=int, default=12,
              help='The maximal depth for paths describing context candidates in the assembly graph')
@click.option('--maximal-gap-ratio', type=float, default=1.5,
              help="The maximal ratio between the length of the gene and the gap between it's contexts in the database")
@click.option('--max-context-len', type=int, default=2500, help='The maximal length for context candidates')
@click.option('--gene-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for locating a gene in the graph')
@click.option('--paths-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for matching a context candidate to a reference sequence')
@click.option('--skip-assembly', is_flag=True, default=False,
              help='A flag that indicates whether or not to skip the assembly step. If the flag is set to True, the argument --assembly--dir must be supplied and direct to the results of a SPAdes run')
def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads, kraken_output_path,
                   reads_ratio_th, metadata_path, references_dir, merged_filtered_fasta,
                   genes_path, depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th,
                   paths_pident_filtering_th, skip_assembly):
    """GInGeR - A tool for analyzing the genomic contexts of genes in metagenomic samples.

    \b
    SHORT_READS_1 - R1 fastq or fastq.gzip file

    \b
    SHORT_READS_2 - R2 fastq or fastq.gzip file

    \b
    GENES_PATH - A fasta file with the genes of interest

    """
    # # filter reference database using kraken
    # if merged_filtered_fasta is None:
    #     merged_filtered_fasta = f'{references_dir}/merged_filtered_ref_db.fasta'
    #     if kraken_output_path is None:
    #         kraken_output_path = f'{out_dir}/kraken_output_file.tsv'
    #     rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads,
    #                                          kraken_output_path, reads_ratio_th,
    #                                          metadata_path, references_dir,
    #                                          merged_filtered_fasta)
    # indexed_reference = sau.generate_index(merged_filtered_fasta, preset=sau.INDEXING_PRESET,
    #                                        just_print=sau.JUST_PRINT_DEFAULT)
    # # run assembly
    # if not skip_assembly:
    #     hau.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_dir, threads)
    # # run tool
    #
    # assembly_graph, genes_with_location_in_graph, records_dict = lg.locate_genes_in_graph(assembly_dir,
    #                                                                                       gene_pident_filtering_th,
    #                                                                                       genes_path, threads, out_dir)
    # if genes_with_location_in_graph is None:
    #     log.info('Genes of interest were not detected in assembly. No results will be generated')
    #     return
    #
    # # get in and out paths
    # in_paths_fasta = c.IN_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    # out_paths_fasta = c.OUT_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    # _ = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, records_dict,
    #                                                           genes_with_location_in_graph, depth_limit,
    #                                                           max_context_len, in_paths_fasta,
    #                                                           out_paths_fasta)
    # map them to the reference
    in_contexts_to_bugs = c.IN_MAPPING_TO_BUGS_PATH_TEMPLATE.format(temp_folder=out_dir)
    out_contexts_to_bugs = c.OUT_MAPPING_TO_BUGS_PATH_TEMPLATE.format(temp_folder=out_dir)
    #
    # sau.map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, indexed_reference, in_contexts_to_bugs,
    #                                    out_contexts_to_bugs, threads)

    # merge and get results
    genes_lengths = extract_genes_lengths(genes_path)
    results = vcc.process_in_and_out_paths_to_results(in_contexts_to_bugs, out_contexts_to_bugs, genes_lengths,
                                                      paths_pident_filtering_th, 0, maximal_gap_ratio)
    context_level_output_path = c.CONTEXT_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    species_level_output_path = c.SPECIES_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    pu.write_context_level_output_to_csv(results, context_level_output_path, metadata_path)
    pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                            species_level_output_path)
    log.info(
        f"GInGer's run completed successfully!\nContext-level output can be found here {context_level_output_path}, species-level output can be found here {species_level_output_path}")


if __name__ == "__main__":
    run_ginger_e2e()
#     # parser = vars(get_command_parser())
#     # arguments that are simply parsed
#     short_reads_1 = parser['short_reads_1']
#     short_reads_2 = parser['short_reads_2']
#     genes_path = parser['genes_path']
#     output_folder = parser['output_folder']
#     long_reads = parser['long_reads']
#     threads = parser['threads']
#     reads_ratio_th = parser['reads_ratio_th']
#     indexed_reference = parser['indexed_reference']
#     depth_limit = parser['depth_limit']
#     maximal_gap_ratio = parser['maximal_gap_ratio']
#     max_context_len = parser['max_context_len']
#     gene_pident_filtering_th = parser['gene_pident_filtering_th']
#     paths_pident_filtering_th = parser['paths_pident_filtering_th']
#     skip_database_filtering = parser['skip_database_filtering']
#     skip_reference_indexing = parser['skip_reference_indexing']
#     skip_assembly = parser['skip_assembly']
#
#     # arguments with more complex parsing
#     assembly_folder = parser['assembly_folder'] if parser['assembly_folder'] is not None else f'{output_folder}/spades'
#     kraken_output_path = f'{output_folder}/ref_db_filtering_kraken_results'
#     references_folder = parser['references_folder'] if parser[
#                                                            'references_folder'] is not None else f'references_folder'
#     metadata_path = parser['metadata_path'] if parser[
#                                                    'metadata_path'] is not None else f'reference_genomes_metadata.tsv'
#     merged_filtered_fasta = f'{output_folder}/merged_filtered_reference_genomes'
#
#     # run ginger
#     ginger_results = run_ginger_e2e(long_reads, short_reads_1, short_reads_2, output_folder, assembly_folder, threads,
#                                     kraken_output_path,
#                                     reads_ratio_th, metadata_path, references_folder, indexed_reference,
#                                     merged_filtered_fasta,
#                                     genes_path,
#                                     depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th,
#                                     paths_pident_filtering_th,
#                                     skip_database_filtering, skip_reference_indexing, skip_assembly)
#
#     context_level_output_path = f'{output_folder}/context_level_matches.csv'
#     species_level_output_path = f'{output_folder}/species_level_matches.csv'
#     pu.write_context_level_output_to_csv(ginger_results, context_level_output_path)
#     pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
#                                                                             species_level_output_path)
#     # with open(f'{output_folder}/output_pickled.pkl', 'wb') as pickle_out_file:
#     #     pickle.dump(ginger_results, pickle_out_file)
