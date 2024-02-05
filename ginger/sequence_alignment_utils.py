from subprocess import run
from pafpy import PafFile
from ginger import pipeline_utils as pu
import os
import datetime as dt
import logging

log = logging.getLogger(__name__)

JUST_PRINT_DEFAULT = False
# MMSEQS:
MMSEQ2_OUTPUT_FORMAT = "'target,query,tstart,tend,nident,qlen'"
MMSEQ2_COMMAND = f"mmseqs easy-search {{query}} {{target}} {{out_file}} mmseqs2_tmp --search-type 2 -a --format-mode 4 --format-output {MMSEQ2_OUTPUT_FORMAT} -c 0.8 --cov-mode 2 --threads {{nthreads}} --mask 0"
MMSEQS_GENES_TO_CONTIGS_HEADER_CONVERSION = {'target': 'contig',
                                             'query': 'gene',
                                             'tstart': 'contig_start',
                                             'tend': 'contig_end',
                                             'nident': 'residue_matches',
                                             'qlen': 'gene_length'}

# TODO do I need to constantly log minimap's output (see assembly_utils) or is it enough to just log it in the end
# minimap:
MINIMAP2_INDEXING_COMMAND = 'minimap2 -x {preset} -d {index_file} {fasta_file}'
MINIMAP2_COMMAND = 'minimap2 -cx {preset} -t {nthreads} {target} {query} > {out_file} -P'
CONTIGS_TO_BUGS_PRESET = 'asm20'
GENES_TO_CONTIGS_PRESET = 'asm20'
ARGS_TO_REFERENCE_CONTIGS = 'asm20'
INDEXING_PRESET = 'asm20'
N_THREADS_DEFAULT = 4

MINIMAP_GENES_TO_CONTIGS_HEADER_CONVERSION = {'query_name': 'gene',
                                              'query_length': 'gene_length',
                                              'query_start': 'gene_start',
                                              'query_end': 'gene_end',
                                              'strand': 'strand',
                                              'target_name': 'contig',
                                              'target_length': 'contig_length',
                                              'target_start': 'contig_start',
                                              'target_end': 'contig_end',
                                              'residue_matches': 'residue_matches',
                                              'alignment_block_length': 'alignment_block_length',
                                              'mapping_quality': 'mapping_quality',
                                              'cg': 'cigar'}

MINIMAP_CONTIGS_TO_REFERENCE_HEADER_CONVERSION = {'query_name': 'contig',
                                                  'query_length': 'contig_length',
                                                  'query_start': 'contig_start',
                                                  'query_end': 'contig_end',
                                                  'strand': 'strand',
                                                  'target_name': 'bug',
                                                  'target_length': 'bug_length',
                                                  'target_start': 'bug_start',
                                                  'target_end': 'bug_end',
                                                  'residue_matches': 'residue_matches',
                                                  'alignment_block_length': 'alignment_block_length',
                                                  'mapping_quality': 'mapping_quality',
                                                  'cg': 'cigar'}

MINIMAP_GT_ARGS_TO_BUGS_HEADER_CONVERSION = {'query_name': 'gene',
                                             'query_length': 'gene_length',
                                             'query_start': 'gene_start',
                                             'query_end': 'gene_end',
                                             'strand': 'strand',
                                             'target_name': 'bug',
                                             'target_length': 'bug_length',
                                             'target_start': 'bug_start',
                                             'target_end': 'bug_end',
                                             'residue_matches': 'residue_matches',
                                             'alignment_block_length': 'alignment_block_length',
                                             'mapping_quality': 'mapping_quality',
                                             }


def generate_index(fasta_file, preset=INDEXING_PRESET, just_print=JUST_PRINT_DEFAULT):
    if fasta_file.endswith('fasta.gz'):
        index_file = fasta_file[:-len('fasta.gz')] + 'mmi'
    elif fasta_file.endswith('fasta.gzip'):
        index_file = fasta_file[:-len('fasta.gzip')] + 'mmi'
    elif fasta_file.endswith('fasta'):
        index_file = fasta_file[:-len('fasta')] + 'mmi'
    else:
        raise Exception('fasta file name should end with fasta.gzip,  fasta.gz or fasta')
    command = MINIMAP2_INDEXING_COMMAND.format(preset=preset, index_file=index_file, fasta_file=fasta_file)
    log.info(f'{dt.datetime.now()} running Minimap2 indexing: {command}')
    if just_print:
        return index_file
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'Minimap2 stderr: {command_output.stderr}')
        raise Exception('Minimap2 failed to run. GInGeR will abort.')
    else:
        logging.info(f'Minimap2 run successfully. stdout: {command_output.stdout}')
    return index_file


def map_contexts_to_bugs(contigs_path, reference_path, contigs_to_bugs_path,
                         just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_minimap2_paf(contigs_path, reference_path, contigs_to_bugs_path, CONTIGS_TO_BUGS_PRESET, just_print,
                             nthreads)


def map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path,
                         just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_mmseqs2(genes_path, contigs_path, genes_to_contigs_path, nthreads)


def map_genes_to_reference(args_path, reference_path, genes_to_reference_path,
                           just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_mmseqs2(args_path, reference_path, genes_to_reference_path, nthreads)


def _run_mmseqs2(query, target, out_file, nthreads=1):
    command = MMSEQ2_COMMAND.format(target=target, query=query, out_file=out_file, nthreads=nthreads)
    logging.info(f'running mmseq2: {command}')
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'mmseq2 stderr: {command_output.stderr}')
    else:
        logging.info(f'mmseq2 run successfully')
    return out_file


def _run_minimap2_paf(query, target, out_file, preset='asm20', just_print=JUST_PRINT_DEFAULT, nthreads=1):
    command = MINIMAP2_COMMAND.format(preset=preset, target=target, query=query, out_file=out_file,
                                      nthreads=nthreads)
    pu.check_and_makedir(out_file)
    logging.info(f'running minimap2: {command}')
    if just_print:
        return out_file
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'minimap2 stderr: {command_output.stderr}')
    else:
        logging.info(f'minimap2 run successfully. stdout: {command_output.stdout}')
    # TODO make this method return the parsed file (had issues with running it on the cluster because I couldn't generate a conda environment)
    return out_file


#
# def filter_genes_in_contigs(df, query_field, pident_filtering_th, alignment_length_filtering_th):
#     # TODO What's the best way to calculate a score for the match?
#     df[f'matches_{query_field}_length_ratio'] = df['residue_matches'] / df[f'{query_field}_length']
#     return df[(df[f'matches_{query_field}_length_ratio'] > pident_filtering_th)]
#
#
# def filter_contigs_in_bugs(df, query_field, pident_filtering_th, alignment_length_filtering_th):
#     df[f'matches_{query_field}_length_ratio'] = df['residue_matches'] / df[
#         f'alignment_block_length']  # note that the name of the column here isn't really correct but I keep it like this to make it easier to generate the bed file
#     return df[(df.alignment_block_length > alignment_length_filtering_th) & (
#             df[f'matches_{query_field}_length_ratio'] > pident_filtering_th)]


def read_and_filter_mmseq2_matches(match_object_constructor: callable, alignment_path: str, pident_filtering_th: float):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        next(f) # skip header
        mmseq2_results = [match_object_constructor(line) for line in f]
    if log.level == logging.DEBUG:
        log.debug(
            f'{dt.datetime.now()} {len(mmseq2_results)} alignments for {len(set([match.gene for match in mmseq2_results]))} genes were found in {alignment_path}')
    filtered_mmseq2_results = [paf_line for paf_line in mmseq2_results if paf_line.match_score > pident_filtering_th]

    if log.level == logging.DEBUG:
        log.debug(
            f'{dt.datetime.now()} {len(filtered_mmseq2_results)} alignments for {len(set([match.gene for match in filtered_mmseq2_results]))} genes were left after filtering ')
    return filtered_mmseq2_results


def read_and_filter_minimap_matches(match_object_constructor: callable, alignment_path: str,
                                    pident_filtering_th: float):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        minimap_results = [match_object_constructor(gene_to_contig) for gene_to_contig in PafFile(f)]

    if log.level == logging.DEBUG:
        minimap_results = minimap_results
        log.debug(
            f'{dt.datetime.now()} {len(minimap_results)} alignments for {len(set([match.gene for match in minimap_results]))} genes were found in {alignment_path}')
    filtered_minimap_results = [paf_line for paf_line in minimap_results if paf_line.match_score > pident_filtering_th]

    if log.level == logging.DEBUG:
        log.debug(
            f'{dt.datetime.now()} {len(filtered_minimap_results)} alignments for {len(set([match.gene for match in filtered_minimap_results]))} genes were left after filtering ')
    return filtered_minimap_results


@pu.step_timing
def map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, reference_path, in_mapping_to_bugs_path,
                                   out_mapping_to_bugs_path, n_minimap_threads):
    map_contexts_to_bugs(in_paths_fasta, reference_path, in_mapping_to_bugs_path, nthreads=n_minimap_threads)
    map_contexts_to_bugs(out_paths_fasta, reference_path, out_mapping_to_bugs_path, nthreads=n_minimap_threads)
