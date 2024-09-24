from pafpy import PafFile
from ginger import pipeline_utils as pu
import os
import datetime as dt
import logging
from tqdm import tqdm
from subprocess import run, Popen, PIPE

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
CONTIGS_TO_REF_GENOMES_PRESET = 'asm20'
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
                                                  'target_name': 'ref_genome',
                                                  'target_length': 'ref_genome_length',
                                                  'target_start': 'ref_genome_start',
                                                  'target_end': 'ref_genome_end',
                                                  'residue_matches': 'residue_matches',
                                                  'alignment_block_length': 'alignment_block_length',
                                                  'mapping_quality': 'mapping_quality',
                                                  'cg': 'cigar'}

MINIMAP_GT_ARGS_TO_BUGS_HEADER_CONVERSION = {'query_name': 'gene',
                                             'query_length': 'gene_length',
                                             'query_start': 'gene_start',
                                             'query_end': 'gene_end',
                                             'strand': 'strand',
                                             'target_name': 'ref_genome',
                                             'target_length': 'ref_genome_length',
                                             'target_start': 'ref_genome_start',
                                             'target_end': 'ref_genome_end',
                                             'residue_matches': 'residue_matches',
                                             'alignment_block_length': 'alignment_block_length',
                                             'mapping_quality': 'mapping_quality',
                                             }


@pu.step_timing
def generate_index(fasta_file, preset):
    # TODO add tqdm
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

    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'Minimap2 stderr: {command_output.stderr}')
        raise Exception('Minimap2 failed to run. GInGeR will abort.')
    else:
        logging.info(f'Minimap2 run successfully. stdout: {command_output.stdout}')
    return index_file


def map_contexts_to_ref_genomes(contigs_path, reference_path, contigs_to_ref_genomes_path,
                                just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_minimap2_paf(contigs_path, reference_path, contigs_to_ref_genomes_path, CONTIGS_TO_REF_GENOMES_PRESET, just_print,
                             nthreads)

def map_nodes_to_contigs(nodes_path, contigs_path, nodes_to_contigs_path,
                         just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_minimap2_paf(nodes_path, contigs_path, nodes_to_contigs_path, CONTIGS_TO_REF_GENOMES_PRESET, just_print,
                             nthreads)


def map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path,
                         just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_mmseqs2(genes_path, contigs_path, genes_to_contigs_path, nthreads)


def map_genes_to_reference(args_path, reference_path, genes_to_reference_path,
                           just_print=JUST_PRINT_DEFAULT, nthreads=N_THREADS_DEFAULT):
    return _run_mmseqs2(args_path, reference_path, genes_to_reference_path, nthreads)


@pu.step_timing
def _run_mmseqs2(query, target, out_file, nthreads=1):
    command = MMSEQ2_COMMAND.format(target=target, query=query, out_file=out_file, nthreads=nthreads)
    pu.check_and_makedir(out_file)
    logging.info(f'running mmseq2: {command}')
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'mmseq2 stderr: {command_output.stderr}')
        raise Exception('mmseq2 failed to run. GInGeR will abort.')
    else:
        logging.info(f'mmseq2 ran successfully')
    return out_file


@pu.step_timing
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
        raise Exception('minimap2 failed to run. GInGeR will abort.')
    else:
        logging.info(f'minimap2 ran successfully')
    return out_file
    # with Popen(command.split(' '), stdout=PIPE) as minimap_process:
    #     output_lines = [output_line for output_line in tqdm(iter(lambda: minimap_process.stdout.readline(), b""))]
    #     if minimap_process.returncode:
    #         logging.error(f'minimap2 stderr: {minimap_process.stderr}')
    #     else:
    #         logging.info(f'minimap2 run successfully. stdout: {minimap_process.stdout}')
    #     return out_file


@pu.step_timing
def read_and_filter_mmseq2_matches(match_object_constructor: callable, alignment_path: str, pident_filtering_th: float):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        next(f)  # skip header
        mmseq2_results = (match_object_constructor(line) for line in f)
        if log.level == logging.DEBUG:  # keeping the generator will be more memory efficient, so convert this to list only if needed for debugging
            mmseq2_results = list(mmseq2_results)
            log.debug(
                f'{dt.datetime.now()} {len(mmseq2_results)} alignments for {len(set([match.gene for match in mmseq2_results]))} genes were found in {alignment_path}')
        filtered_mmseq2_results = [paf_line for paf_line in mmseq2_results if paf_line.score > pident_filtering_th]
        if log.level == logging.DEBUG:
            log.debug(
                f'{dt.datetime.now()} {len(filtered_mmseq2_results)} alignments for {len(set([match.gene for match in filtered_mmseq2_results]))} genes were left after filtering ')
    return filtered_mmseq2_results

@pu.step_timing
def read_and_filter_minimap_matches(match_object_constructor: callable, alignment_path: str,
                                    pident_filtering_th: float, *argv):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        minimap_results = [match_object_constructor(gene_to_contig, *argv) for gene_to_contig in PafFile(f)]

    if log.level == logging.DEBUG:
        minimap_results = minimap_results
        log.debug(
            f'{dt.datetime.now()} {len(minimap_results)} alignments for {len(set([match.gene for match in minimap_results]))} genes were found in {alignment_path}')
    filtered_minimap_results = [paf_line for paf_line in minimap_results if paf_line.score > pident_filtering_th]

    if log.level == logging.DEBUG:
        log.debug(
            f'{dt.datetime.now()} {len(filtered_minimap_results)} alignments for {len(set([match.gene for match in filtered_minimap_results]))} genes were left after filtering ')
    return filtered_minimap_results


@pu.step_timing
def map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, reference_path, in_mapping_to_ref_genomes_path,
                                   out_mapping_to_ref_genomes_path, n_minimap_threads):
    map_contexts_to_ref_genomes(in_paths_fasta, reference_path, in_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
    map_contexts_to_ref_genomes(out_paths_fasta, reference_path, out_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
