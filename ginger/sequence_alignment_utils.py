from pafpy import PafFile
from ginger import pipeline_utils as pu
import os
import datetime as dt
import logging
from tqdm import tqdm
from subprocess import run, Popen, PIPE
import itertools

log = logging.getLogger(__name__)

IOU_TH = 0.3

JUST_PRINT_DEFAULT = False
# MMSEQS:
MMSEQ2_OUTPUT_FORMAT = "'target,query,tstart,tend,nident,qlen'"
MMSEQ2_COMMAND = f"mmseqs easy-search {{query}} {{target}} {{out_file}} {{temp_dir}} --search-type 2 -a --format-mode 4 --format-output {MMSEQ2_OUTPUT_FORMAT} -c 0.8 --cov-mode 2 --threads {{nthreads}} --mask 0"
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
    log.info(f'Running Minimap2 indexing - {command}')

    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'Minimap2 indexing failed: {command_output.stderr}')
        raise Exception('Minimap2 indexing failed - GInGeR aborted')
    else:
        logging.info(f'Minimap2 indexing completed successfully')
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

@pu.step_timing
def _run_mmseqs2(query, target, out_file, nthreads=1):
    # temp dir should be under the dir of the out file
    command = MMSEQ2_COMMAND.format(target=target, query=query, out_file=out_file, temp_dir=f'{os.path.dirname(out_file)}/mmseqs_tmp', nthreads=nthreads)
    pu.check_and_makedir(out_file)
    logging.info(f'Running MMseqs2 - {command}')
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'MMseqs2 failed: {command_output.stderr}')
        raise Exception('MMseqs2 failed - GInGeR aborted')
    else:
        logging.info(f'MMseqs2 completed successfully')
    return out_file


@pu.step_timing
def _run_minimap2_paf(query, target, out_file, preset='asm20', just_print=JUST_PRINT_DEFAULT, nthreads=1):
    command = MINIMAP2_COMMAND.format(preset=preset, target=target, query=query, out_file=out_file,
                                      nthreads=nthreads)
    pu.check_and_makedir(out_file)
    logging.info(f'Running Minimap2 - {command}')
    if just_print:
        return out_file
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        logging.error(f'Minimap2 failed: {command_output.stderr}')
        raise Exception('Minimap2 failed - GInGeR aborted')
    else:
        logging.info(f'Minimap2 completed successfully')
    return out_file
    # with Popen(command.split(' '), stdout=PIPE) as minimap_process:
    #     output_lines = [output_line for output_line in tqdm(iter(lambda: minimap_process.stdout.readline(), b""))]
    #     if minimap_process.returncode:
    #         logging.error(f'minimap2 stderr: {minimap_process.stderr}')
    #     else:
    #         logging.info(f'minimap2 run successfully. stdout: {minimap_process.stdout}')
    #     return out_file
def interval_iou(start1, end1, start2, end2):
    try:
        intersection = max(0, min(end1, end2) - max(start1, start2))
        if intersection == 0:
            return 0
        union = max(end1, end2) - min(start1, start2)
        output = intersection / union
    except Exception as e:
        log.error(f'{e} - start1: {start1}, end1: {end1}, start2: {start2}, end2: {end2}')
        raise e
    return output


def is_similar_to_representatives(representatives, gene_paths_to_ref_genome_match, iou_th):
    for rep in representatives:
        if interval_iou(gene_paths_to_ref_genome_match.start, gene_paths_to_ref_genome_match.end, rep.start,
                        rep.end) > iou_th:
            return True
    return False

def non_max_suppresion(matches_dict):
    nms_dict = {}
    for bug, matches in matches_dict.items():
        nmsd_list = non_max_suppression_single_class(matches)
        if nmsd_list:
            nms_dict[bug] = nmsd_list
    return nms_dict


def non_max_suppression_single_class(matches, sorting_func=lambda x: (x.score, x.end - x.start, x.start, x.gene),
                                     iou_th=IOU_TH,
                                     score_th=0):
    # to maintain consistencythe soring function should put the best matches first, and if there are matches with the same score, sort them by their start position, and if they are the same, sort them by the gene name.
    sorted_matches = sorted(matches, key=sorting_func, reverse=True)
    representative_matches = []
    for match in sorted_matches:
        if not is_similar_to_representatives(representative_matches, match, iou_th):
            representative_matches.append(match)
    return representative_matches

@pu.step_timing
def read_and_filter_mmseq2_matches(match_object_constructor: callable, alignment_path: str, pident_filtering_th: float, nms=True):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        next(f)  # skip header
        mmseq2_results = [match_object_constructor(line) for line in f]
        if log.level == logging.DEBUG:  # keeping the generator will be more memory efficient, so convert this to list only if needed for debugging
            mmseq2_results = list(mmseq2_results)
            log.debug(
                f'Found {len(mmseq2_results)} alignments for {len(set([match.gene for match in mmseq2_results]))} genes in {alignment_path}')
        
        # filter and group by contig
        filtered_mmseqs_by_contig_dict = {}
        for match in mmseq2_results:
            if match.score > pident_filtering_th:
                if match.contig in filtered_mmseqs_by_contig_dict:
                    filtered_mmseqs_by_contig_dict[match.contig].append(match)
                else:
                    filtered_mmseqs_by_contig_dict[match.contig] = [match]
        
        if log.level == logging.DEBUG:
            log.debug(
                f'Retained {len([matches for matches in filtered_mmseqs_by_contig_dict.values()])} alignments for species and {len(set([match.gene for matches in filtered_mmseqs_by_contig_dict.values() for match in matches]))} genes after filtering')

        if nms:
            filtered_mmseqs_by_contig_dict = non_max_suppresion(filtered_mmseqs_by_contig_dict)
            
            if log.level == logging.DEBUG:
                log.debug(
                    f'Retained {len([matches for matches in filtered_mmseqs_by_contig_dict.values()])} alignments for species and {len(set([match.gene for matches in filtered_mmseqs_by_contig_dict.values() for match in matches]))} genes after taking to gene per lofaction')
        
        # Flatten dictionary values into a single list
        filtered_matches_list = list(itertools.chain.from_iterable(filtered_mmseqs_by_contig_dict.values()))
    return filtered_matches_list

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
            f'Found {len(minimap_results)} alignments for {len(set([match.gene for match in minimap_results]))} genes in {alignment_path}')
    filtered_minimap_results = [paf_line for paf_line in minimap_results if paf_line.score > pident_filtering_th]

    if log.level == logging.DEBUG:
        log.debug(
            f'Retained {len(filtered_minimap_results)} alignments for {len(set([match.gene for match in filtered_minimap_results]))} genes after filtering')
    return filtered_minimap_results


@pu.step_timing
def map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, reference_path, in_mapping_to_ref_genomes_path,
                                   out_mapping_to_ref_genomes_path, n_minimap_threads):
    map_contexts_to_ref_genomes(in_paths_fasta, reference_path, in_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
    map_contexts_to_ref_genomes(out_paths_fasta, reference_path, out_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
