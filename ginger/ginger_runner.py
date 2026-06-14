import sys
import logging
import click
from Bio import SeqIO
import os
import pandas as pd
import shutil
from glob import glob
import pickle

from ginger import locating_genes_in_graph as lg
from ginger import reference_database_utils as rdu
from ginger import assembly_utils as au
from ginger import extract_contexts_candidates as ecc
from ginger import sequence_alignment_utils as sau
from ginger import verify_context_candidates as vcc
from ginger import pipeline_utils as pu
from ginger import constants as c


logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

def cleanup_intermediate_files(out_dir, keep_options):
    """Remove intermediate files based on keep_options."""
    if 'all' in keep_options:
        return
    
    log.info('Cleaning up intermediate files')
    
    # Define file patterns for each category
    cleanup_map = {
        'assembly': ['SPAdes'],
        'alignment': ['*.paf', '*.m8', 'mmseqs_tmp', 'nodes_to_contigs_w_gaps.paf'],
        'sequences': ['all_in_paths.fasta', 'all_out_paths.fasta'],
        'kraken': ['kraken_*.tsv', 'bracken_*.tsv'],
        'reference': ['merged_filtered_ref_db.*', 'references_used.csv']
    }
    
    # Remove categories not in keep_options
    # 'final' is not a real category, it just means "keep final CSV results"
    categories_to_remove = [cat for cat in cleanup_map.keys() if cat not in keep_options]
    
    for category in categories_to_remove:
        patterns = cleanup_map[category]
        for pattern in patterns:
            # Handle both files and directories
            if '*' in pattern:
                # Glob pattern
                matches = glob(os.path.join(out_dir, pattern))
                for match in matches:
                    try:
                        if os.path.isdir(match):
                            shutil.rmtree(match)
                            log.debug(f'Removed directory {match}')
                        else:
                            os.remove(match)
                            log.debug(f'Removed file {match}')
                    except Exception as e:
                        log.warning(f'Failed to remove {match}: {e}')
            else:
                # Exact file/directory name
                path = os.path.join(out_dir, pattern)
                if os.path.exists(path):
                    try:
                        if os.path.isdir(path):
                            shutil.rmtree(path)
                            log.debug(f'Removed directory {path}')
                        else:
                            os.remove(path)
                            log.debug(f'Removed file {path}')
                    except Exception as e:
                        log.warning(f'Failed to remove {path}: {e}')


@click.command()
@click.argument('short-reads-1', required=True, type=click.Path(exists=True))
@click.argument('short-reads-2', required=True, type=click.Path(exists=True))
@click.argument('genes-path', required=True, type=click.Path(exists=True))
@click.argument('out-dir', required=True, type=click.Path())
@click.option('--long-reads', type=click.Path(exists=True), help='A fastq or fastq.gzip file of Oxford Nanopore reads')
@click.option('--assembly-dir', default=None,
              help="Specifies where to save the assembly results. In case of pre-ran assembly, please insert the path do the spades output directory")
@click.option('--threads', '-t', type=int, default=1,
              help='Number of threads that will be used for running Kraken2, SPAdes and Minimap2')
@click.option('--kraken-output-path', default=None, help="A path for saving Kraken2's output")
@click.option('--kraken-db', type=click.Path(),
              default=os.path.join(os.path.dirname(__file__), '..', 'kraken2_db_uhgg_v2.0.2'),
              help='The path to UHGG\'s Kraken2 database directory')
@click.option('--species-coverage-threshold', type=float, default=10,
              help='The minimal estimated sequencing coverage required for including a species in the analysis. Coverage is estimated as: bracken_estimated_reads * (avg_len_R1 + avg_len_R2) / median_genome_length, where median genome length is computed from the top references per species (by Quality) capped by --max-species-representatives. Default 10.')
@click.option('--max-species-representatives', type=int, default=100,
              help='The maximal references per species that will be downloaded from UHGG and taken into account in the aggregation of results at the species level')
@click.option('--reference-genomes-metadata', type=click.Path(),
              default=os.path.join(os.path.dirname(__file__), 'UHGG-metadata.tsv'),
              help='The path to the reference database metadata table')
@click.option('--downloaded-references-dir', type=click.Path(), default='references_dir',
              help='The directory to which GInGeR will download missing reference genomes from UHGG. This folder can be shared for all runs of GInGer in order to avoid the same file being  downloaded and saved multiple times')
@click.option('--sample-specific-references', type=click.Path(), default=None,
              help='A fasta, fasta.gz or mmi (minimap indexed) file that will be used a reference database (using this will skip the stages of creating a sample specific database based on the species detected in the sample by Kraken2)')
@click.option('--depth-limit', type=int, default=12,
              help='The maximal depth for paths describing context candidates in the assembly graph')
@click.option('--max-gap-ratio', type=float, default=1.5,
              help="The maximal ratio between the length of the gene and the gap between it's contexts in the database")
@click.option('--max-context-len', type=int, default=2500, help='The maximal length for context candidates')
@click.option('--min-context-len', type=int, default=0, help='The minimal length for context candidates')
@click.option('--gene-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for locating a gene in the graph')
@click.option('--paths-pident-filtering-th', type=float, default=0.9,
              help='The minimal % of matched base pairs required for matching a context candidate to a reference sequence')
@click.option('--keep-intermediate', multiple=True,
              type=click.Choice(['all', 'final', 'assembly', 'alignment', 'sequences', 'kraken', 'reference'],
                               case_sensitive=False),
              default=['all'],
              help='Specify which intermediate files to keep. Options: all (default, keep everything), final (only result CSVs), assembly (SPAdes output), alignment (PAF/M8 files), sequences (FASTA files), kraken (Kraken2/Bracken output), reference (reference database files). Can specify multiple by repeating the flag: --keep-intermediate final --keep-intermediate assembly')
@click.option('--skip-assembly', is_flag=True, default=False,
              help='A flag that indicates whether or not to skip the assembly step. If the flag is set to True, the argument --assembly--dir must be supplied and direct to the results of a SPAdes run')
@click.option('--return-all-gene-matches', is_flag=True, default=False,
              help='By default, GInGeR applies non-max-suppression (NMS) to the alignment of genes of interest to the assembly graph, keeping only top scoring matches and removing redundant overlapping matches. If this flag is set to True, all gene matches will be returned without applying NMS.')
@click.option('--nms-iou-threshold', type=float, default=0.8,
              help='The IoU (Intersection over Union) threshold used for non-max-suppression when filtering overlapping gene matches. Gene matches with IoU > this threshold are considered overlapping. Only used when --return-all-gene-matches is False. Default: 0.8')
def run_ginger_e2e(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads, kraken_output_path,
                   kraken_db, species_coverage_threshold, reference_genomes_metadata, downloaded_references_dir, sample_specific_references, genes_path, depth_limit,
                   max_gap_ratio, max_context_len, min_context_len, gene_pident_filtering_th,
                   paths_pident_filtering_th, keep_intermediate, skip_assembly, max_species_representatives, return_all_gene_matches, nms_iou_threshold):
    """GInGeR - A tool for analyzing the genomic contexts of genes in metagenomic samples.

    \b
    SHORT_READS_1 - R1 fastq or fastq.gzip file

    \b
    SHORT_READS_2 - R2 fastq or fastq.gzip file
t
    \b
    GENES_PATH - A fasta file with the genes of interest

    \b
    OUT_DIR - A path specifying where to save GInGeR's output

    """
    return ginger_e2e_func(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads, kraken_output_path,
                           kraken_db, species_coverage_threshold, reference_genomes_metadata, downloaded_references_dir, sample_specific_references, genes_path,
                           depth_limit, max_gap_ratio, min_context_len, max_context_len, gene_pident_filtering_th,
                           paths_pident_filtering_th, keep_intermediate, skip_assembly, max_species_representatives, return_all_gene_matches, nms_iou_threshold)


def ginger_e2e_func(long_reads, short_reads_1, short_reads_2, out_dir, assembly_dir, threads, kraken_output_path,
                    kraken_db, species_coverage_threshold, reference_genomes_metadata, downloaded_references_dir, sample_specific_references, genes_path, depth_limit,
                    max_gap_ratio, min_context_len, max_context_len, gene_pident_filtering_th,
                    paths_pident_filtering_th, keep_intermediate, skip_assembly, max_species_representatives, return_all_gene_matches, nms_iou_threshold):
    # Log the command that was run
    log.info(f"Running GInGeR with command: {' '.join(sys.argv)}")
    
    # create output directory if it doesn't exist
    pu.check_and_make_dir_no_file_name(out_dir)
    # filter reference database using kraken
    references_used_path = f'{out_dir}/references_used.csv'
    if sample_specific_references is None:
        sample_specific_references = f'{out_dir}/merged_filtered_ref_db.fasta'
        # if the file was not specified or the specified file does not exist
        if kraken_output_path is None or not os.path.exists(kraken_output_path):
            kraken_output_path = f'{out_dir}/kraken_output_file.tsv'
            kraken_report_path = f'{out_dir}/kraken_report_file.tsv'
            bracken_output = f'{out_dir}/bracken_output_file.tsv'
            bracken_report = f'{out_dir}/bracken_report_file.tsv'
            rdu.get_filtered_references_database(short_reads_1, short_reads_2, threads, kraken_output_path,
                                                 kraken_report_path, bracken_output, bracken_report, species_coverage_threshold,
                                                 reference_genomes_metadata, downloaded_references_dir, sample_specific_references,
                                                 references_used_path,
                                                 max_species_representatives, kraken_db)
    if not sample_specific_references.endswith('mmi'):
        indexed_reference = sau.generate_index(sample_specific_references, sau.INDEXING_PRESET)
    else:
        indexed_reference = sample_specific_references
    # run assembly
    if assembly_dir is None:
        assembly_dir = f'{out_dir}/SPAdes'
    if not skip_assembly:
        au.run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, assembly_dir, threads)
    # run tool

    assembly_graph, genes_with_location_in_graph, assembly_graph_nodes = lg.locate_genes_in_graph(assembly_dir,
                                                                                                  gene_pident_filtering_th,
                                                                                                  genes_path,
                                                                                                  threads,
                                                                                                  out_dir,
                                                                                                  return_all_gene_matches,
                                                                                                  nms_iou_threshold)
    if not genes_with_location_in_graph:
        log.info(
            'No genes of interest detected in assembly. GInGeR run stopped - no results generated')
        return

    genes_detected_no_species_match_output_path = c.GENES_DETECTED_IN_GRAPH_WITH_NO_SPECIES_MATCH_OUTPUT_TEMPLATE.format(
        out_dir=out_dir
    )

    # get in and out paths
    in_paths_fasta = c.IN_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    out_paths_fasta = c.OUT_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
    gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, assembly_graph_nodes,
                                                              genes_with_location_in_graph, depth_limit,
                                                              min_context_len, max_context_len, in_paths_fasta,
                                                              out_paths_fasta)

    # map them to the reference
    in_contexts_to_ref_genomes = c.IN_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE.format(temp_folder=out_dir)
    out_contexts_to_ref_genomes = c.OUT_MAPPING_TO_REF_GENOMES_PATH_TEMPLATE.format(temp_folder=out_dir)

    sau.map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, indexed_reference, in_contexts_to_ref_genomes,
                                       out_contexts_to_ref_genomes, threads)

    # merge and get results
    context_level_results = vcc.process_in_and_out_paths_to_results(in_contexts_to_ref_genomes,
                                                                    out_contexts_to_ref_genomes,
                                                                    gene_lengths, paths_pident_filtering_th, 0,
                                                                    max_gap_ratio, reference_genomes_metadata)
    if not context_level_results:
        pu.write_genes_detected_in_graph_with_no_species_match(
            genes_with_location_in_graph,
            matched_genes=set(),
            csv_path=genes_detected_no_species_match_output_path,
        )
        log.info(
            'No matching pairs of incoming and outgoing contexts found in reference sequences. GInGeR run stopped - no results generated')
        return
    context_level_output_path = c.CONTEXT_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    species_level_output_path = c.SPECIES_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    subspecies_level_output_path = c.SUBSPECIES_LEVEL_OUTPUT_TEMPLATE.format(out_dir=out_dir)
    pu.write_context_level_output_to_csv(context_level_results, context_level_output_path, reference_genomes_metadata)
    species_level_df = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path,
                                                                                               reference_genomes_metadata,
                                                                                               species_level_output_path,
                                                                                               max_species_representatives)

    if os.path.exists(references_used_path) and 'subspecies' in pd.read_table(references_used_path).columns:
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path,
                                                                                references_used_path,
                                                                                subspecies_level_output_path,
                                                                                max_species_representatives,
                                                                                'subspecies')

    matched_genes = set()
    if species_level_df is not None and len(species_level_df.index) > 0:
        # species_level_df has a MultiIndex (gene, species)
        matched_genes = set(species_level_df.index.get_level_values(0))

    pu.write_genes_detected_in_graph_with_no_species_match(
        genes_with_location_in_graph,
        matched_genes=matched_genes,
        csv_path=genes_detected_no_species_match_output_path,
    )
    # Clean up intermediate files if requested
    cleanup_intermediate_files(out_dir, keep_intermediate)
    
    log.info(
        f"GInGeR completed successfully. Context-level output: {context_level_output_path}, Species-level output: {species_level_output_path}")

if __name__ == "__main__":
    run_ginger_e2e()
