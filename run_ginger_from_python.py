import ginger
import pickle

if __name__ == "__main__":
    LAB ='/sci/labs/morani/morani/icore-data/lab'
    PROJECT_DIR = f'{LAB}/Projects/hybrid_assembler'
    PIPELINE_TREE = f'{PROJECT_DIR}/pipeline_tree'

    long_reads = f'{PIPELINE_TREE}/zymo_standard/data/long_reads/ERR3152364_GridION_sequencing_sample_250000.fastq.gz'
    short_reads_1 = f'{PIPELINE_TREE}/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_1_sample_1000000.fastq.gz'
    short_reads_2 = f'{PIPELINE_TREE}/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_2_sample_1000000.fastq.gz'
    output_folder = f'{LAB}/Projects/hybrid_assembler/test_run'
    assembly_folder = None
    threads = 16
    kraken_output_path = f'{output_folder}/kraken_folder/kraken_output_file.tsv'
    reads_ratio_th = 0.01
    metadata_path = f'{LAB}/Tools/UnifiedHumanGastrointestinalGenome/genomes-all_metadata.tsv'
    references_folder = f'{output_folder}/references_folder/'
    indexed_reference = None
    merged_filtered_fasta = f'{output_folder}/kraken_folder/merged_filtered_ref_db.fasta'
    genes_path = f'{PROJECT_DIR}/pipeline_tree/genes_from_amir_erez/processed/genes_without_dup_locations.fasta'
    depth_limit = 12
    maximal_gap_ratio = 1.5
    max_context_len = 2500
    gene_pident_filtering_th = 0.9
    paths_pident_filtering_th = 0.9
    skip_database_filtering = False
    skip_reference_indexing = False
    skip_assembly = False


    ginger_results = ginger.run_ginger_e2e(long_reads, short_reads_1, short_reads_2, output_folder, assembly_folder,
                                           threads,
                                           kraken_output_path,
                                           reads_ratio_th, metadata_path, references_folder, indexed_reference,
                                           merged_filtered_fasta,
                                           genes_path,
                                           depth_limit, maximal_gap_ratio, max_context_len, gene_pident_filtering_th,
                                           paths_pident_filtering_th,
                                           skip_database_filtering, skip_reference_indexing, skip_assembly)
    with open(f'{output_folder}/output_pickled.pkl', 'wb') as pickle_out_file:
        pickle.dump(ginger_results, pickle_out_file)
