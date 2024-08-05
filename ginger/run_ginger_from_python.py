from ginger import ginger_runner
import pickle
import shutil

if __name__ == "__main__":
    params = dict(
        short_reads_1='/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R1.fastq.gz',
        # '/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R1_1M_sample.fastq',
        short_reads_2='/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R2.fastq.gz',
        # '/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads/MSA-1006_1_R2_1M_sample.fastq',
        # genes_path='/sci/labs/morani/morani/icore-data/lab/Data/BacterialGenes/fn_MH0327_GL0091324.fasta',
        # '/sci/labs/morani/morani/icore-data/lab/Data/BacterialGenes/IGC_proteins_nms_on_atcc_1006.fasta',
        # '/sci/labs/morani/morani/icore-data/lab/Data/BacterialGenes/IGC_proteins_no_dups.fasta',
        genes_path= '/sci/labs/morani/morani/icore-data/lab/Data/BacterialGenes/sample_IGC_proteins_no_dups.fasta',
        # '/sci/labs/morani/morani/icore-data/lab/Data/BacterialGenes/sample_IGC_proteins_no_dups_pos_strain.fasta',
        out_dir='ginger_out_dir',
        long_reads='/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/nanopore/SRR9847864.fastq.gz',
        assembly_dir='/sci/backup/morani/lab/Projects/GInGeR/analysis_data/SPAdes',
        threads=2,
        kraken_output_path=None,
        reads_ratio_th=0.01,
        max_species_representatives=100,
        metadata_path='/sci/backup/morani/lab/Projects/GInGeR/GInGeR/ginger/UHGG-metadata.tsv',
        references_dir='/sci/backup/morani/lab/Projects/GInGeR/analysis_data/analysis_references_dir',
        merged_filtered_fasta='/sci/backup/morani/lab/Projects/GInGeR/analysis_data/analysis_references_dir/MGYG000005036.fasta',
        # merged_filtered_fasta='/sci/backup/morani/lab/Projects/GInGeR/analysis_data/analysis_references_dir/merged_filtered_ref_db.mmi',
        depth_limit=12,
        maximal_gap_ratio=1.5,
        max_context_len=2500,
        gene_pident_filtering_th=0.9,
        paths_pident_filtering_th=0.9,
        skip_assembly=True)

    # remove out_dir if exists
    shutil.rmtree(params['out_dir'], ignore_errors=True)
    ginger_results = ginger_runner.ginger_e2e_func(**params)
    # with open(f'{out_dir}/output_pickled.pkl', 'wb') as pickle_out_file:
    #     pickle.dump(ginger_results, pickle_out_file)
