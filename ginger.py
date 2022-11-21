import sys
import logging
import click

import locating_genes_in_graph as lg
import reference_database_utils as rdu
import hybrid_assembly_utils as hau
import extract_contexts_candidates as ecc
import sequence_alignment_utils as sau
import verify_context_candidates as vcc
import constants as c

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


# def extract_genes_lengths(genes_path):
#     gene_lengths = {}
#     with open(genes_path) as f:
#         fasta_sequences = SeqIO.parse(f, 'fasta')
#         for s in fasta_sequences:
#             gene_length = len(s.seq)
#             gene_name = s.id
#
#             gene_lengths[gene_name] = gene_length
#     return gene_lengths
#

# kraken_output_path = f'{output_folder}/kraken_folder/kraken_output_file.tsv'
# reads_ratio_th = 0.01
# metadata_path = f'{LAB}/Tools/UnifiedHumanGastrointestinalGenome/genomes-all_metadata.tsv'
# references_folder = f'{output_folder}/references_folder/'
# indexed_reference = None
# merged_filtered_fasta = f'{output_folder}/kraken_folder/merged_filtered_ref_db.fasta'

@click.command()
@click.option('short_reads_1', required=True, type=click.Path(exists=True), help='R1 fastq or fastq.gzip file')
@click.option('short_reads_2', required=True, type=click.Path(exists=True), help='R2 fastq or fastq.gzip file')
@click.option('genes_path', required=True, type=click.Path(exists=True), help='A fasta file with the genes of interest')
@click.option('long_reads', type=click.Path(exists=True), help='A fastq or fastq.gzip file of Oxford Nanopore reads')
@click.option('-o', '--out-folder', default='ginger_output', help="A path specifying where to save GInGeR's output",
              type=click.Path())
@click.option('assembly_dir', default='ginger_output/SPAdes',
              help="Specifies where to save the assembly results. In case of pre-ran assembly, please insert the path do the spades output directory",
              type=click.Path())
@click.option('threads', type=int, default=1,
              help='Number of threads that will be used for running Kraken2, SPAdes and Minimap2')
@click.option('kraken_output_path', default='ginger_output/kraken_folder/kraken_output_file.tsv', type=click.Path())
@click.option('reads_ratio_th', type=float, default=0.9)
@click.option('metadata_path')
@click.option('references_folder')
@click.option('indexed_reference')
@click.option('merged_filtered_fasta')
@click.option('depth_limit')
@click.option('maximal_gap_ratio')
@click.option('max_context_len')
@click.option('gene_pident_filtering_th')
@click.option('paths_pident_filtering_th')
@click.option('skip_database_preprocessing')
@click.option('skip_assembly')
def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, out_folder, assembly_dir, threads, kraken_output_path,
                   reads_ratio_th, metadata_path, references_folder, indexed_reference, merged_filtered_fasta,
                   genes_path, depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th,
                   paths_pident_filtering_th,
                   skip_database_preprocessing, skip_assembly):
    # filter reference database using kraken
    if not skip_database_preprocessing:
        rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path, reads_ratio_th,
                                             metadata_path, references_folder, merged_filtered_fasta)
        indexed_reference = sau.generate_index(merged_filtered_fasta, preset=sau.INDEXING_PRESET,
                                               just_print=sau.JUST_PRINT_DEFAULT)
    # run assembly
    if not skip_assembly:
        if assembly_dir is None:
            assembly_dir = f'{out_folder}/{c.ASSEMBLY_FOLDER_NAME}'
        hau.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_dir, threads)
    # run tool

    assembly_graph, genes_to_contigs, records_dict = lg.locate_genes_in_graph(assembly_dir, gene_pident_filtering_th,
                                                                              genes_path, threads, out_folder)

    # TODO if there are no in paths, there is no need to calc out paths. but maybe we should always have in paths. check this
    # get in and out paths
    in_paths_fasta = c.IN_PATHS_FASTA_TEMPLATE.format(temp_folder=out_folder)
    out_paths_fasta = c.OUT_PATHS_FASTA_TEMPLATE.format(temp_folder=out_folder)
    genes_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, records_dict,
                                                                          genes_to_contigs, depth_limit,
                                                                          max_context_len, in_paths_fasta,
                                                                          out_paths_fasta)
    # map them to the reference
    in_contexts_to_bugs = c.IN_MAPPING_TO_BUGS_PATH_TEMPLATE.format(temp_folder=out_folder)
    out_contexts_to_bugs = c.OUT_MAPPING_TO_BUGS_PATH_TEMPLATE.format(temp_folder=out_folder)

    sau.map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, indexed_reference, in_contexts_to_bugs,
                                       out_contexts_to_bugs, threads)

    # merge and get results
    # genes_lengths = extract_genes_lengths(genes_path)
    results = vcc.process_in_and_out_paths_to_results(in_contexts_to_bugs, out_contexts_to_bugs, genes_lengths,
                                                      paths_pident_filtering_th, 0, maximal_gap_ratio)
    return results

# if __name__ == "__main__":
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
