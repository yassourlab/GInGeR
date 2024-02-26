from ginger import ginger_runner
import pickle
 # TODO this file should be deleted
if __name__ == "__main__":
    LAB = '/sci/labs/morani/morani/icore-data/lab'
    # LAB = '/Volumes/netta.barak/lab'
    PROJECT_DIR = f'{LAB}/Projects/hybrid_assembler'
    PIPELINE_TREE = f'{PROJECT_DIR}/pipeline_tree'

    long_reads = f'{PIPELINE_TREE}/zymo_standard/data/long_reads/ERR3152364_GridION_sequencing_sample_250000.fastq.gz'
    short_reads_1 = f'{PIPELINE_TREE}/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_1_sample_1000000.fastq.gz'
    short_reads_2 = f'{PIPELINE_TREE}/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_2_sample_1000000.fastq.gz'
    genes_path = f'{PROJECT_DIR}/pipeline_tree/genes_from_amir_erez/processed/genes_without_dup_locations.fasta'
    out_dir = f'{LAB}/Projects/hybrid_assembler/test_run'
    assembly_dir = None
    threads = 16
    kraken_output_path = f'{out_dir}/kraken_folder/kraken_output_file.tsv'
    reads_ratio_th = 0.01
    metadata_path = f'{LAB}/Tools/UnifiedHumanGastrointestinalGenome/genomes-all_metadata.tsv'
    references_folder = f'{out_dir}/references_folder/'
    indexed_reference = None
    merged_filtered_fasta = f'{out_dir}/kraken_folder/merged_filtered_ref_db.fasta'

    depth_limit = 12
    maximal_gap_ratio = 1.5
    max_context_len = 2500
    gene_pident_filtering_th = 0.9
    paths_pident_filtering_th = 0.9
    skip_database_filtering = False
    skip_reference_indexing = False
    skip_assembly = False

    ginger_results = ginger_runner.run_ginger_e2e(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads,
                                               kraken_output_path,
                                               reads_ratio_th, metadata_path, references_folder, merged_filtered_fasta,
                                               genes_path, depth_limit, maximal_gap_ratio, max_context_len,
                                               gene_pident_filtering_th,
                                               paths_pident_filtering_th, skip_assembly)
    with open(f'{out_dir}/output_pickled.pkl', 'wb') as pickle_out_file:
        pickle.dump(ginger_results, pickle_out_file)
