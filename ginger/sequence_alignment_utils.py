from pafpy import PafFile
from ginger import pipeline_utils as pu
import os
import logging
from collections import defaultdict
import itertools

log = logging.getLogger(__name__)

IOU_TH = 0.3

# MMSEQS:
MMSEQ2_OUTPUT_FORMAT = "'target,query,tstart,tend,nident,qlen'"
MMSEQ2_COMMAND = f"mmseqs easy-search {{query}} {{target}} {{out_file}} {{temp_dir}} --search-type 2 -a --format-mode 4 --format-output {MMSEQ2_OUTPUT_FORMAT} -c 0.8 --cov-mode 2 --threads {{nthreads}} --mask 0"

# minimap:
MINIMAP2_INDEXING_COMMAND = 'minimap2 -x {preset} -d {index_file} {fasta_file}'
MINIMAP2_COMMAND = 'minimap2 -cx {preset} -t {nthreads} {target} {query} > {out_file} -P'
CONTIGS_TO_REF_GENOMES_PRESET = 'asm20'
INDEXING_PRESET = 'asm20'
N_THREADS_DEFAULT = 4
# 'fasta' has to come last - it is a suffix of the other two
FASTA_SUFFIXES = ('fasta.gz', 'fasta.gzip', 'fasta')


@pu.step_timing
def generate_index(fasta_file, preset):
    for suffix in FASTA_SUFFIXES:
        if fasta_file.endswith(suffix):
            index_file = fasta_file[:-len(suffix)] + 'mmi'
            break
    else:
        raise Exception('fasta file name should end with fasta.gzip,  fasta.gz or fasta')
    pu.run_tool('Minimap2 indexing',
                MINIMAP2_INDEXING_COMMAND.format(preset=preset, index_file=index_file, fasta_file=fasta_file))
    return index_file


@pu.step_timing
def run_minimap2_paf(query, target, out_file, nthreads=N_THREADS_DEFAULT, preset=CONTIGS_TO_REF_GENOMES_PRESET):
    pu.check_and_makedir(out_file)
    pu.run_tool('Minimap2', MINIMAP2_COMMAND.format(preset=preset, target=target, query=query, out_file=out_file,
                                                    nthreads=nthreads))
    return out_file


@pu.step_timing
def map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path, nthreads=N_THREADS_DEFAULT):
    # temp dir should be under the dir of the out file
    command = MMSEQ2_COMMAND.format(query=genes_path, target=contigs_path, out_file=genes_to_contigs_path,
                                    temp_dir=f'{os.path.dirname(genes_to_contigs_path)}/mmseqs_tmp', nthreads=nthreads)
    pu.check_and_makedir(genes_to_contigs_path)
    pu.run_tool('MMseqs2', command)
    return genes_to_contigs_path


def interval_iou(start1, end1, start2, end2):
    intersection = max(0, min(end1, end2) - max(start1, start2))
    if intersection == 0:
        return 0
    return intersection / (max(end1, end2) - min(start1, start2))


def is_similar_to_representatives(representatives, gene_paths_to_ref_genome_match, iou_th):
    for rep in representatives:
        if interval_iou(gene_paths_to_ref_genome_match.start, gene_paths_to_ref_genome_match.end, rep.start,
                        rep.end) > iou_th:
            return True
    return False


def non_max_suppression_single_class(matches, sorting_func=lambda x: (x.score, x.end - x.start, x.gene),
                                     iou_th=IOU_TH, reverse=True):
    """Keeps the best of every group of matches that overlap each other by more than iou_th.

    sorting_func and reverse together have to put the best match first, and break ties the match length and the gene name.
    """
    sorted_matches = sorted(matches, key=sorting_func, reverse=reverse)
    representative_matches = []
    for match in sorted_matches:
        if not is_similar_to_representatives(representative_matches, match, iou_th):
            representative_matches.append(match)
    return representative_matches


@pu.step_timing
def read_and_filter_mmseq2_matches(match_object_constructor: callable, alignment_path: str, pident_filtering_th: float,
                                   nms=True, nms_iou_threshold=IOU_TH):
    if os.path.getsize(alignment_path) == 0:
        return None
    matches_by_contig = defaultdict(list)
    with open(alignment_path, 'r') as f:
        next(f)  # skip header
        for line in f:
            match = match_object_constructor(line)
            if match.score > pident_filtering_th:
                matches_by_contig[match.contig].append(match)

    if nms:
        # per contig, so that two copies of a gene on different contigs never suppress each other
        matches_by_contig = {contig: non_max_suppression_single_class(matches, iou_th=nms_iou_threshold)
                             for contig, matches in matches_by_contig.items()}
    filtered_matches_list = list(itertools.chain.from_iterable(matches_by_contig.values()))
    log.debug(f'kept {len(filtered_matches_list)} alignments for '
              f'{len({match.gene for match in filtered_matches_list})} genes from {alignment_path}')
    return filtered_matches_list


@pu.step_timing
def read_and_filter_minimap_matches(match_object_constructor: callable, alignment_path: str,
                                    pident_filtering_th: float, *match_constructor_args):
    if os.path.getsize(alignment_path) == 0:
        return None
    with open(alignment_path, 'r') as f:
        matches = [match_object_constructor(paf_line, *match_constructor_args) for paf_line in PafFile(f)]

    filtered_minimap_results = [match for match in matches if match.score > pident_filtering_th]
    log.debug(f'kept {len(filtered_minimap_results)} of {len(matches)} alignments for '
              f'{len({match.gene for match in filtered_minimap_results})} genes from {alignment_path}')
    return filtered_minimap_results


@pu.step_timing
def map_in_and_out_contexts_to_ref(in_paths_fasta, out_paths_fasta, reference_path, in_mapping_to_ref_genomes_path,
                                   out_mapping_to_ref_genomes_path, n_minimap_threads):
    run_minimap2_paf(in_paths_fasta, reference_path, in_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
    run_minimap2_paf(out_paths_fasta, reference_path, out_mapping_to_ref_genomes_path, nthreads=n_minimap_threads)
