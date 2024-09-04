from ginger import ginger_runner
import pickle
import shutil

if __name__ == "__main__":
    params = {
        'short_reads_1': '/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R1_1M_sample.fastq',
        'short_reads_2': '/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R2_1M_sample.fastq',
        'genes_path': '/sci/backup/morani/lab/Projects/GInGeR/analysis_data/250twins_37276_GL0178535.fasta',
        'out_dir': 'ginger_out_dir', 'long_reads': None,
        'assembly_dir': '/sci/backup/morani/lab/Projects/GInGeR/analysis_data/SPAdes_short_3M_long_0',
        'kraken_output_path': None, 'reads_ratio_th': 0.01, 'max_species_representatives': 200,
        'metadata_path': '/sci/labs/morani/morani/icore-data/lab/Data/UnifiedHumanGastrointestinalGenome/genomes-partial_metadata_with_species.tsv',
        'references_dir': '/sci/backup/morani/lab/Projects/GInGeR/analysis_data/analysis_references_dir_200',
        'merged_filtered_fasta': '/sci/backup/morani/lab/Projects/GInGeR/analysis_data/analysis_references_dir_200/merged_filtered_ref_db_add_atcc.mmi',
        'depth_limit': 12, 'maximal_gap_ratio': 1.5, 'max_context_len': 2500, 'min_context_len': 2000,
        'gene_pident_filtering_th': 0.9, 'paths_pident_filtering_th': 0.9, 'skip_assembly': True, 'threads': '2'}

    # remove out_dir if exists
    shutil.rmtree(params['out_dir'], ignore_errors=True)
    ginger_results = ginger_runner.ginger_e2e_func(**params)
    # with open(f'{out_dir}/output_pickled.pkl', 'wb') as pickle_out_file:
    #     pickle.dump(ginger_results, pickle_out_file)
