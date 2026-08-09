import pyfastg
from collections import Counter
from ginger import pipeline_utils as pu
import logging
from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import constants as c
from Bio import SeqIO
from typing import Dict, Iterator, List, Set
import pandas as pd
import os

log = logging.getLogger(__name__)

def get_node_without_adj(long_node_name):
    split_by_dots = long_node_name.split(':')[0]
    split_by_comma_dot = long_node_name.split(';')[0]
    if len(split_by_comma_dot) <= len(split_by_dots):
        return split_by_comma_dot
    return split_by_dots


def get_short_node_name(long_node_name):
    node_without_adj = get_node_without_adj(long_node_name)
    node_num = node_without_adj.split('_')[1]
    last_char_chuku = node_without_adj[-1] == "'"
    if last_char_chuku:
        return node_num + '-'
    else:
        return node_num + '+'


def add_location_in_graph_based_on_contigs_paths(nodes_sequences_dict, nodes_in_path, gene_contig_match, offset=0):
    #  this is not the exact start but it's good enough
    #  offset is where nodes_in_path starts in contig coordinates - 0 for a contig assembled from a
    #  single graph path, the segment's origin for one of the segments of a gap-containing contig
    contig_start = gene_contig_match.start - offset
    contig_end = gene_contig_match.end - offset
    if len(nodes_in_path) == 1:
        gene_contig_match.nodes_list = nodes_in_path
        gene_contig_match.start_in_first_node = contig_start
    else:
        prev_seq = str(nodes_sequences_dict[nodes_in_path[0]].seq)
        end = len(prev_seq)
        if contig_start <= end:  # in case that the match starts in the first node
            start_in_first_node = contig_start
            nodes_for_genes = [nodes_in_path[0]]
        else:
            nodes_for_genes = []
            start_in_first_node = None
        for node in nodes_in_path[1:]:
            cur_seq = str(nodes_sequences_dict[node].seq)
            try:
                k = pu.get_sequence_overlap(prev_seq, cur_seq)
            except Exception as e:
                log.error(f'error! {str(e)} {node} ')
                raise Exception
            start = end - k
            end = start + len(cur_seq)
            if pu.intervals_overlap(start, end, contig_start, contig_end):
                if start_in_first_node is None and start <= contig_start <= end:
                    start_in_first_node = contig_start - start
                nodes_for_genes.append(node)
            else:
                if nodes_for_genes:
                    break
            prev_seq = cur_seq
        gene_contig_match.nodes_list = nodes_for_genes
        gene_contig_match.start_in_first_node = start_in_first_node
    # TODO I'm modifying and then returning the same object. I think it's not the best practice
    return gene_contig_match


def node_oriented_with_contig(node_name, strand):
    """The short name of the graph node that runs in the same direction as the contig.

    A nodes_list is always oriented with the contig - that is what contigs.paths gives the routes
    that read it - so that the contexts extracted from it flank the gene the way they flank it on
    the contig. An alignment can reach that node from either of its two fastg records: the forward
    one aligning on the plus strand, or the reverse complement one aligning on the minus strand.
    Both orientations of every edge are nodes of the graph, so flipping is always possible.
    """
    short_node_name = get_short_node_name(node_name)
    if strand == '+':
        return short_node_name
    return short_node_name[:-1] + ('-' if short_node_name.endswith('+') else '+')


def get_start_in_first_node(gene_contig_match, top_nodes_to_contigs_match):
    """Where the gene starts inside the node. The node runs with the contig, so its coordinates are
    the contig's, shifted to where it aligned."""
    return gene_contig_match.start - top_nodes_to_contigs_match.contig_start


def add_location_in_graph_based_on_nodes_to_contigs(gene_contig_match, nodes_to_contigs_df):
    """Locates a gene on the single graph node that aligned over the gene's start, for when the
    contig's path could not place it.

    The node has to cover the start rather than merely overlap the gene, because start_in_first_node
    is an offset into it.
    """
    covers_gene_start = ((nodes_to_contigs_df.contig_start <= gene_contig_match.start) &
                         (gene_contig_match.start <= nodes_to_contigs_df.contig_end))
    filtered_nodes_to_contigs = nodes_to_contigs_df[
        (nodes_to_contigs_df.contig == gene_contig_match.contig) & covers_gene_start]

    if filtered_nodes_to_contigs.empty:
        return None
    top_nodes_to_contigs_match = filtered_nodes_to_contigs.iloc[0]
    gene_contig_match.nodes_list = [node_oriented_with_contig(top_nodes_to_contigs_match.node,
                                                              top_nodes_to_contigs_match.strand)]
    gene_contig_match.start_in_first_node = get_start_in_first_node(gene_contig_match, top_nodes_to_contigs_match)
    return gene_contig_match


def anchor_segment_in_contig(segment, contig_nodes_to_contigs):
    """Where a path segment of a gap-containing contig starts in contig coordinates, taken from the
    alignment of the segment's first node to the contig. Gap lengths are unknown, so this can't be
    derived from node lengths. Returns None if the first node has no alignment above the score cutoff
    applied in map_nodes_to_contigs_w_gaps.
    """
    # a node and its reverse complement align to the same contig interval and only one of them
    # survives the deduplication in map_nodes_to_contigs_w_gaps, so the node's orientation in the
    # segment is ignored here and only the node number is matched
    first_node_number = segment[0][:-1]
    anchors = contig_nodes_to_contigs[
        contig_nodes_to_contigs.node.apply(lambda node: node.split('_')[1]) == first_node_number]
    if anchors.empty:
        return None
    return anchors.iloc[0].contig_start


def add_location_in_graph_based_on_gappy_contig_paths(nodes_sequences_dict, segments, gene_contig_match,
                                                      nodes_to_contigs_df):
    """Locates a gene in the graph when its contig was assembled from several graph paths joined
    using paired-end evidence, by walking the single segment that covers the gene. Returns None if
    no segment covering the gene could be anchored in contig coordinates.
    """
    contig_nodes_to_contigs = nodes_to_contigs_df[nodes_to_contigs_df.contig == gene_contig_match.contig]
    for segment in segments:
        origin = anchor_segment_in_contig(segment, contig_nodes_to_contigs)
        if origin is None:
            continue
        segment_length = len(pu.generate_str_from_list_of_nodes(nodes_sequences_dict, segment)[0])
        if origin <= gene_contig_match.start and gene_contig_match.end <= origin + segment_length:
            return add_location_in_graph_based_on_contigs_paths(nodes_sequences_dict, segment, gene_contig_match,
                                                                offset=origin)
    return None


def add_node_list_to_genes_to_contigs(genes_to_contigs: Iterator[mc.GeneContigMatch],
                                      parsed_paths: Dict[str, List[List[str]]],
                                      nodes_sequences_dict: Dict[str, SeqIO.SeqRecord],
                                      nodes_to_contigs_df: pd.DataFrame,
                                      contigs_with_gaps: Set[str] = frozenset()):
    """Adds a nodes_list and a start_in_first_node to every gene-contig match that can be located in
    the assembly graph, and returns the matches worth analyzing.

    A match on one of contigs_with_gaps is kept even when it could not be located, because its
    context is taken from the contig rather than from the graph.
    """
    matches_to_analyze = []
    genes_located_per_route = Counter()

    for gene_contig_match in genes_to_contigs:
        segments = parsed_paths.get(gene_contig_match.contig, [])
        located_gene_contig_match = None
        route = None
        if len(segments) == 1:
            located_gene_contig_match = add_location_in_graph_based_on_contigs_paths(nodes_sequences_dict,
                                                                                     segments[0],
                                                                                     gene_contig_match)
            route = 'the contig path'
        elif segments:
            located_gene_contig_match = add_location_in_graph_based_on_gappy_contig_paths(nodes_sequences_dict,
                                                                                          segments,
                                                                                          gene_contig_match,
                                                                                          nodes_to_contigs_df)
            route = 'a gap-containing contig path segment'
        if located_gene_contig_match is None:
            located_gene_contig_match = add_location_in_graph_based_on_nodes_to_contigs(gene_contig_match,
                                                                                        nodes_to_contigs_df)
            route = 'a single aligned node'
        if gene_contig_match.start_in_first_node is not None:
            genes_located_per_route[route] += 1
        if located_gene_contig_match is not None or gene_contig_match.contig in contigs_with_gaps:
            matches_to_analyze.append(gene_contig_match)
    for route, count in genes_located_per_route.most_common():
        log.info(f'located {count} genes in the assembly graph by {route}')
    return matches_to_analyze


def get_nodes_dict_from_fastg_file(assembly_graph_path: str) -> Dict[str, SeqIO.SeqRecord]:
    """
    parses an assembly graph fastg file (reads it as a fasta file)
    :param assembly_graph_path: a fastg file representing the assembly graph (one of the outputs of spades)
    :return:
    """
    with open(assembly_graph_path) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    nodes_sequences_dict = {get_short_node_name(record.name): record for record in records}
    return nodes_sequences_dict


@pu.step_timing
def find_genes_in_contigs(genes_path: str, contigs_path: str, n_minimap_threads: int,
                          pident_filtering_th: float, genes_to_contigs_path: str,
                          return_all_gene_matches: bool = False, nms_iou_threshold: float = 0.8) -> Iterator[mc.GeneContigMatch]:
    if not os.path.exists(genes_to_contigs_path) or not os.path.isfile(genes_to_contigs_path):
        log.info(f'running mmseqs2 to find genes in contigs')
        genes_to_contigs_path = sau.map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path,
                                                         nthreads=n_minimap_threads)
    genes_to_contigs = sau.read_and_filter_mmseq2_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                          pident_filtering_th, nms=not return_all_gene_matches,
                                                          nms_iou_threshold=nms_iou_threshold)
    
    return genes_to_contigs



@pu.step_timing
def locate_genes_in_graph(assembly_dir: str, gene_pident_filtering_th: float, genes_path: str, n_threads: int,
                          temp_folder: str, return_all_gene_matches: bool = False, nms_iou_threshold: float = 0.8):  # -> Tuple[networkx.DiGraph,??? ,Dict[str, SeqIO.SeqRecord]]
    contigs_path = c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    assembly_graph_path = c.ASSEMBLY_GRAPH_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    genes_to_contigs_path = c.GENES_TO_CONTIGS_TEMPLATE.format(temp_files_path=temp_folder)
    nodes_to_contigs_w_gaps_path = c.NODES_TO_CONTIGS_W_GAPS_TEMPLATE.format(temp_files_path=temp_folder)

    # find genes in contigs
    genes_to_contigs = find_genes_in_contigs(genes_path, contigs_path, n_threads,
                                             gene_pident_filtering_th, genes_to_contigs_path,
                                             return_all_gene_matches, nms_iou_threshold)
    if not genes_to_contigs:
        return None, None, None, None

    assembly_graph_nodes = get_nodes_dict_from_fastg_file(assembly_graph_path)
    assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    # find nodes in contigs with gaps
    parsed_paths, contigs_with_gaps = pu.parse_paths_file(
        c.PATHS_PATH_TEMPLATE.format(assembly_dir=assembly_dir), assembly_graph.nodes)
    nodes_to_contigs_df = map_nodes_to_contigs_w_gaps(contigs_with_gaps, assembly_graph_path, contigs_path, n_threads,
                                                      nodes_to_contigs_w_gaps_path)

    genes_to_analyze = add_node_list_to_genes_to_contigs(genes_to_contigs, parsed_paths, assembly_graph_nodes,
                                                         nodes_to_contigs_df, contigs_with_gaps)
    n_located = sum(1 for match in genes_to_analyze if match.start_in_first_node is not None)
    log.info(f'found locations in the assembly graph for {n_located} genes, and kept '
             f'{len(genes_to_analyze) - n_located} more genes that were found on gap-containing contigs')
    return assembly_graph, genes_to_analyze, assembly_graph_nodes, contigs_with_gaps


@pu.step_timing
def map_nodes_to_contigs_w_gaps(contigs_with_gaps, assembly_graph_path, contigs_path, n_threads,
                                nodes_to_contigs_w_gaps_path):
    # filter contigs fasta to keep only contigs with gaps
    contigs_w_gaps_path = contigs_path.replace('.fasta', '_w_gaps.fasta')
    SeqIO.write([contig for contig in SeqIO.parse(contigs_path, 'fasta') if contig.id in contigs_with_gaps],
                contigs_w_gaps_path, 'fasta')
    # find nodes in contigs with gaps
    nodes_to_contigs_path = sau.map_nodes_to_contigs(assembly_graph_path, contigs_w_gaps_path,
                                                     nodes_to_contigs_w_gaps_path,
                                                     nthreads=n_threads)
    nodes_to_contigs_df = pu.minimap_results_from_path(nodes_to_contigs_path)
    node_to_contig_columns = ['contig', 'contig_start', 'contig_end', 'node', 'score', 'strand']
    if len(nodes_to_contigs_df):
        nodes_to_contigs_df['node'] = nodes_to_contigs_df.qname.apply(lambda x: x[:-1].split(':')[0])
        nodes_to_contigs_df['score'] = nodes_to_contigs_df.mlen / nodes_to_contigs_df.qlen
        nodes_to_contigs_df = \
            nodes_to_contigs_df[(nodes_to_contigs_df.score > 0.95)].rename(
                columns={'tname': 'contig', 'tstart': 'contig_start', 'tend': 'contig_end'})[
                node_to_contig_columns]
        nodes_to_contigs_df.sort_values(['score', 'node'], inplace=True, ascending=(False, True))
        nodes_to_contigs_df.drop_duplicates(subset=['contig', 'contig_start', 'contig_end', 'score'],
                                            keep='first', inplace=True)
    else:  # keep the columns so that downstream filtering works on an empty result too
        nodes_to_contigs_df = pd.DataFrame(columns=node_to_contig_columns)
    return nodes_to_contigs_df
