CONTIGS_PATH_TEMPLATE = "{assembly_dir}/contigs.fasta"
PATHS_PATH_TEMPLATE = "{assembly_dir}/contigs.paths"
ASSEMBLY_GRAPH_PATH_TEMPLATE = "{assembly_dir}/assembly_graph.fastg"
IN_PATHS_FASTA_TEMPLATE = "{temp_folder}/all_in_paths.fasta"
OUT_PATHS_FASTA_TEMPLATE = "{temp_folder}/all_out_paths.fasta"
IN_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE = '{temp_folder}/in_paths_to_reference.paf'
OUT_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE = '{temp_folder}/out_paths_to_reference.paf'
GENES_TO_CONTIGS_TEMPLATE = '{temp_files_path}/genes_to_contigs.m8'
NODES_TO_CONTIGS_W_GAPS_TEMPLATE = '{temp_files_path}/nodes_to_contigs_w_gaps.paf'
CONTEXT_LEVEL_OUTPUT_TEMPLATE = '{out_dir}/context_level_matches.csv'
SPECIES_LEVEL_OUTPUT_TEMPLATE = '{out_dir}/species_level_matches.csv'
SUBSPECIES_LEVEL_OUTPUT_TEMPLATE = '{out_dir}/subspecies_level_matches.csv'
GENES_DETECTED_IN_GRAPH_WITH_NO_SPECIES_MATCH_OUTPUT_TEMPLATE = '{out_dir}/genes_detected_in_graph_with_no_species_match.csv'
# a result rather than an intermediate - the sequence behind every context level row, for the user to
# take to a genome browser, and GeNomad's input on the way
IN_GENE_OUT_CONTEXTS_FASTA_TEMPLATE = '{out_dir}/in_gene_out_contexts.fasta'
PLASMID_SUMMARY_TEMPLATE = '{out_dir}/plasmid_summary.tsv'
GENOMAD_OUTPUT_DIR_TEMPLATE = '{out_dir}/genomad_output'
SPECIES_INCLUDED_IN_ANALYSIS_TEMPLATE = '{out_dir}/species_included_in_analysis.csv'

