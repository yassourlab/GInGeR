import sys
from Bio import SeqIO
import logging
import argparse
import pandas as pd
import pickle
import locate_genes as lg

import pipeline_utils as pu
import reference_database_utils as rdu
import hybrid_assembly_utils as hau
import process_context_candidates as iopa
import extract_contexts_candidates as ecc
import sequence_alignment_utils as sau
import process_context_candidates as pcc
import constants

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)


#
# def get_command_parser():
#     # TODO - deal with conda env
#     # TODO make sure parsing from terminal works well
#     parser = argparse.ArgumentParser(prog='GInGeR',
#                                      description='A tool for analyzing the genomic regions of genes of interest in a metagenomic sample')
#     # must haves
#     parser.add_argument('short_reads_1')
#     parser.add_argument('short_reads_2')
#     parser.add_argument('genes_path')
#     parser.add_argument('output_folder')
#
#     # TODO make it possible to work with other databases other than the UHGG Kraken one
#     # reference genomes database preprocessing optional parameters:
#     parser.add_argument('--skip-species-downloading', action=argparse.BooleanOptionalAction)
#     parser.add_argument('--skip-reference-indexing', action=argparse.BooleanOptionalAction)
#     parser.add_argument('--minimal-reads-ratio', default=0.01, type=float)
#     parser.add_argument('--reference-metadata-path')
#
#     parser.add_argument('--reference-genome-db')
#
#     # Assembly optional parameters:
#     parser.add_argument('--skip-assembly', action=argparse.BooleanOptionalAction)
#     parser.add_argument('--long-reads')
#     parser.add_argument('--assembly-folder')
#
#     # GInGeR optional parameters:
#     parser.add_argument('--depth-limit', default=12, type=int)
#     parser.add_argument('--maximal-gap-ratio', default=1.5, type=float)
#     parser.add_argument('--max-context-len', default=2500, type=int)
#     parser.add_argument('--gene-pident-filtering-th', default=0.9, type=float)
#     parser.add_argument('--paths-pident-filtering-th', default=0.9, type=float)
#     parser.add_argument('--threads', default=1, type=int)
#     return parser
#
#
#
#
#
# def run_in_out_paths_approach(assembly_dir, genes_path, reference_path,
#                               n_minimap_threads, depth_limit, maximal_gap_ratio, context_len,
#                               gene_pident_filtering_th,
#                               paths_pident_filtering_th, temp_folder):


def extract_genes_lengths(genes_path):
    gene_lengths = {}
    with open(genes_path) as f:
        fasta_sequences = SeqIO.parse(f, 'fasta')
        for s in fasta_sequences:
            gene_length = len(s.seq)
            gene_name = s.id

            gene_lengths[gene_name] = gene_length
    return gene_lengths

def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, temp_folder, assembly_dir, threads, kraken_output_path,
                   reads_ratio_th, metadata_path, references_folder, indexed_reference, merged_filtered_fasta,
                   genes_path,
                   depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th, paths_pident_filtering_th,
                   skip_database_filtering, skip_reference_indexing, skip_assembly):
    # filter reference database using kraken
    if not skip_database_filtering:
        pu.check_and_makedir(kraken_output_path)
        pu.check_and_make_dir_no_file_name(references_folder)
        rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path, reads_ratio_th,
                                             metadata_path, references_folder, merged_filtered_fasta)
    if not skip_reference_indexing:
        indexed_reference = sau.generate_index(merged_filtered_fasta, preset=sau.INDEXING_PRESET,
                                               just_print=sau.JUST_PRINT_DEFAULT)  # TODO can multithread this also

    # run assembly
    if assembly_dir is None:
        assembly_dir = f'{temp_folder}/{constants.ASSEMBLY_FOLDER_NAME}'
    if not skip_assembly:
        hau.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_dir, threads)

    # run tool

    contigs_path = f"{assembly_dir}/{'contigs.fasta'}"
    paths_path = f"{assembly_dir}/{'contigs.paths'}"
    assembly_graph_path = f"{assembly_dir}/{'assembly_graph.fastg'}"
    in_paths_fasta = f"{temp_folder}/all_in_paths.fasta"
    out_paths_fasta = f"{temp_folder}/all_out_paths.fasta"
    in_mapping_to_bugs_path = f'{temp_folder}/in_paths_to_reference.paf'
    out_mapping_to_bugs_path = f'{temp_folder}/out_paths_to_reference.paf'

    assembly_graph, genes_to_contigs, records_dict = lg.find_genes_in_graph(assembly_graph_path, contigs_path,
                                                                            gene_pident_filtering_th, genes_path,
                                                                            threads, paths_path, temp_folder)
    # TODO if there are no in paths, there is no need to calc out paths. but maybe we should always have in paths. check this
    # get in and out paths
    genes_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, records_dict,
                                                                          genes_to_contigs, depth_limit,
                                                                          max_context_len, in_paths_fasta,
                                                                          out_paths_fasta)
    # map them to the reference
    sau.map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, indexed_reference, in_mapping_to_bugs_path,
                                       out_mapping_to_bugs_path, threads)

    # merge and get results
    genes_lengths = extract_genes_lengths(genes_path)
    results = pcc.process_in_and_out_paths_to_results(in_mapping_to_bugs_path, out_mapping_to_bugs_path, genes_lengths,
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
