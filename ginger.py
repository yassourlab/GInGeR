import reference_database_utils as rdu
import hybrid_assembly_utils as hau
import in_out_paths_approach as iopa
import sequence_alignment_utils as sau
import timeit

RUNTIME_PRINTS_PATTERN = '$$$$$$$$$$'
ASSEMBLY_FOLDER_NAME = 'hybrid_spades'


# TODO write parser for all parameters
def define_parser():
    pass
    # Kraken parameters

    # Assembly parameters

    # my tool parameters


# TODO use logging instead of prints
def step_timing(start, step_name):
    stop = timeit.default_timer()
    print(f'{RUNTIME_PRINTS_PATTERN} {step_name} took {(stop - start) / 60} minutes {RUNTIME_PRINTS_PATTERN}')
    return timeit.default_timer()


def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, temp_folder, assembly_folder, threads, kraken_output_path,
                   reads_ratio_th, metadata_path, references_folder, indexed_reference, merged_filtered_fasta,
                   genes_path,
                   depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th, paths_pident_filtering_th,
                   skip_database_filtering, skip_reference_indexing, skip_assemblying):
    start = timeit.default_timer()
    # filter reference database using kraken
    if not skip_database_filtering:
        rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path, reads_ratio_th,
                                             metadata_path, references_folder, merged_filtered_fasta)
    start = step_timing(start, 'database filtering')
    if not skip_reference_indexing:
        indexed_reference = sau.generate_index(merged_filtered_fasta, preset=sau.INDEXING_PRESET,
                                               just_print=sau.JUST_PRINT_DEFAULT)  # TODO can multithread this also

    start = step_timing(start, 'indexing reference database')

    # run assembly
    if assembly_folder is None:
        assembly_folder = f'{temp_folder}/{ASSEMBLY_FOLDER_NAME}'
    if not skip_assemblying:
        hau.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_folder, threads)

    start = step_timing(start, 'running assembly')

    # run tool
    results = iopa.run_in_out_paths_approach(assembly_folder, genes_path, indexed_reference, threads, depth_limit,
                                             maximal_gap_ratio,
                                             max_context_len, gene_pident_filtering_th, paths_pident_filtering_th,
                                             temp_folder)
    step_timing(start, 'finding genes contexts')
    return results
