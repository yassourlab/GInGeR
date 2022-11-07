import sys
import timeit
import logging
import argparse
import pandas as pd
import pickle

import pipeline_utils as pu
import reference_database_utils as rdu
import hybrid_assembly_utils as hau
import process_context_candidates as iopa
import sequence_alignment_utils as sau
import constants

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
log = logging.getLogger(__name__)

RUNTIME_PRINTS_PATTERN = '$$$$$$$$$$'


def get_command_parser():
    # TODO - deal with conda env
    # TODO make sure parsing from terminal works well
    parser = argparse.ArgumentParser(prog='GInGeR',
                                     description='A tool for analyzing the genomic regions of genes of interest in a metagenomic sample')
    # must haves
    parser.add_argument('short_reads_1')
    parser.add_argument('short_reads_2')
    parser.add_argument('genes_path')
    parser.add_argument('output_folder')

    # TODO make it possible to work with other databases other than the UHGG Kraken one
    # reference genomes database preprocessing optional parameters:
    parser.add_argument('--skip-species-downloading', action=argparse.BooleanOptionalAction)
    parser.add_argument('--skip-reference-indexing', action=argparse.BooleanOptionalAction)
    parser.add_argument('--minimal-reads-ratio', default=0.01, type=float)
    parser.add_argument('--reference-metadata-path')

    parser.add_argument('--reference-genome-db')

    # Assembly optional parameters:
    parser.add_argument('--skip-assembly', action=argparse.BooleanOptionalAction)
    parser.add_argument('--long-reads')
    parser.add_argument('--assembly-folder')

    # GInGeR optional parameters:
    parser.add_argument('--depth-limit', default=12, type=int)
    parser.add_argument('--maximal-gap-ratio', default=1.5, type=float)
    parser.add_argument('--max-context-len', default=2500, type=int)
    parser.add_argument('--gene-pident-filtering-th', default=0.9, type=float)
    parser.add_argument('--paths-pident-filtering-th', default=0.9, type=float)
    parser.add_argument('--threads', default=1, type=int)
    return parser


def step_timing(start, step_name):
    stop = timeit.default_timer()
    log.info(f'{RUNTIME_PRINTS_PATTERN} {step_name} took {(stop - start) / 60} minutes {RUNTIME_PRINTS_PATTERN}')
    return timeit.default_timer()


def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, temp_folder, assembly_folder, threads, kraken_output_path,
                   reads_ratio_th, metadata_path, references_folder, indexed_reference, merged_filtered_fasta,
                   genes_path,
                   depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th, paths_pident_filtering_th,
                   skip_database_filtering, skip_reference_indexing, skip_assembly):
    start = timeit.default_timer()
    # filter reference database using kraken
    if not skip_database_filtering:
        pu.check_and_makedir(kraken_output_path)
        pu.check_and_make_dir_no_file_name(references_folder)
        rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path, reads_ratio_th,
                                             metadata_path, references_folder, merged_filtered_fasta)
        start = step_timing(start, 'database filtering')
    if not skip_reference_indexing:
        indexed_reference = sau.generate_index(merged_filtered_fasta, preset=sau.INDEXING_PRESET,
                                               just_print=sau.JUST_PRINT_DEFAULT)  # TODO can multithread this also

        start = step_timing(start, 'indexing reference database')

    # run assembly
    if assembly_folder is None:
        assembly_folder = f'{temp_folder}/{constants.ASSEMBLY_FOLDER_NAME}'
    if not skip_assembly:
        hau.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_folder, threads)

        start = step_timing(start, 'running assembly')

    # run tool
    results = iopa.run_in_out_paths_approach(assembly_folder, genes_path, indexed_reference, threads, depth_limit,
                                             maximal_gap_ratio,
                                             max_context_len, gene_pident_filtering_th, paths_pident_filtering_th,
                                             temp_folder)
    # TODO need to write the output to files
    step_timing(start, 'finding genes contexts')
    return results


def write_context_level_output_to_csv(output, csv_path):
    output_lines = [
        'gene,reference_contig,in_context,out_context,in_context_start,gene_start,gene_end,out_context_end,match_score\n']
    for gene_species_tuple, matches_list in output.items():
        gene, species = gene_species_tuple
        for match in matches_list:
            in_context = match.in_path.query_name
            out_context = match.out_path.query_name
            in_context_start = match.in_path.bug_start
            gene_start = match.start
            gene_end = match.end
            out_context_end = match.out_path.bug_end
            match_score = match.match_score
            output_lines.append(
                f'{gene},{species},{in_context},{out_context},{in_context_start},{gene_start},{gene_end},{out_context_end},{match_score}\n')

    with open(csv_path, 'w') as f:
        f.writelines(output_lines)


def aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                         species_level_output_path):
    context_level_df = pd.read_csv(context_level_output_path)
    context_level_df['Genome'] = context_level_df['reference_contig'].apply(lambda x: x.split('_')[0].split('.')[0])
    metadata_df = pd.read_csv(metadata_path, sep='\t')
    context_level_df_with_metadata = pd.merge(context_level_df, metadata_df[['Genome', 'Species_rep']], on='Genome',
                                              how='left')
    genomes_per_species = metadata_df.Species_rep.value_counts()
    agg_output = context_level_df_with_metadata.groupby(['gene', 'Species_rep']).aggregate(
        {'Genome': ['count'], 'match_score': ['max']})
    agg_output.columns = ['_'.join(col) for col in agg_output.columns.values]
    agg_output = agg_output.merge(genomes_per_species, left_on='Species_rep', right_index=True, how='left')
    agg_output['references_ratio'] = agg_output['Genome_count'] / agg_output['Species_rep']
    agg_output[['references_ratio', 'match_score_max']].to_csv(species_level_output_path)


if __name__ == "__main__":
    parser = vars(get_command_parser())
    # arguments that are simply parsed
    short_reads_1 = parser['short_reads_1']
    short_reads_2 = parser['short_reads_2']
    genes_path = parser['genes_path']
    output_folder = parser['output_folder']
    long_reads = parser['long_reads']
    threads = parser['threads']
    reads_ratio_th = parser['reads_ratio_th']
    indexed_reference = parser['indexed_reference']
    depth_limit = parser['depth_limit']
    maximal_gap_ratio = parser['maximal_gap_ratio']
    max_context_len = parser['max_context_len']
    gene_pident_filtering_th = parser['gene_pident_filtering_th']
    paths_pident_filtering_th = parser['paths_pident_filtering_th']
    skip_database_filtering = parser['skip_database_filtering']
    skip_reference_indexing = parser['skip_reference_indexing']
    skip_assembly = parser['skip_assembly']

    # arguments with more complex parsing
    assembly_folder = parser['assembly_folder'] if parser['assembly_folder'] is not None else f'{output_folder}/spades'
    kraken_output_path = f'{output_folder}/ref_db_filtering_kraken_results'
    references_folder = parser['references_folder'] if parser[
                                                           'references_folder'] is not None else f'references_folder'
    metadata_path = parser['metadata_path'] if parser[
                                                   'metadata_path'] is not None else f'reference_genomes_metadata.tsv'
    merged_filtered_fasta = f'{output_folder}/merged_filtered_reference_genomes'

    # run ginger
    ginger_results = run_ginger_e2e(long_reads, short_reads_1, short_reads_2, output_folder, assembly_folder, threads,
                                    kraken_output_path,
                                    reads_ratio_th, metadata_path, references_folder, indexed_reference,
                                    merged_filtered_fasta,
                                    genes_path,
                                    depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th,
                                    paths_pident_filtering_th,
                                    skip_database_filtering, skip_reference_indexing, skip_assembly)

    context_level_output_path = f'{output_folder}/context_level_matches.csv'
    species_level_output_path = f'{output_folder}/species_level_matches.csv'
    write_context_level_output_to_csv(ginger_results, context_level_output_path)
    aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                         species_level_output_path)
    # with open(f'{output_folder}/output_pickled.pkl', 'wb') as pickle_out_file:
    #     pickle.dump(ginger_results, pickle_out_file)
