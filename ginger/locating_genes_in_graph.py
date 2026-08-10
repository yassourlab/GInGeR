import pyfastg
from collections import Counter, defaultdict
from itertools import groupby
from ginger import pipeline_utils as pu
import logging
from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import constants as c
from Bio import SeqIO
from typing import Dict, Iterator, List, Set
import os

log = logging.getLogger(__name__)

# a node has to be almost entirely covered by its alignment to place it on a contig - a partial hit
# is a repeat or a shared prefix, not the node sitting there
NODE_TO_CONTIG_SCORE_TH = 0.95
# how far apart two alignments may put the same segment before they are treated as disagreeing
# rather than as the same placement seen twice. minimap2 clips alignment ends, so candidates from
# different nodes of one segment differ by tens of bases even when all of them are right
MAX_ANCHOR_DISAGREEMENT = 300


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


def locate_gene_on_path(geometry, nodes_in_path, gene_contig_match, origin=0):
    """Places a gene on a path of the graph, in the path's own coordinates.

    origin is where the path starts in contig coordinates - 0 for a contig assembled from a single
    graph path, the segment's origin for one of the segments of a gap-containing contig.

    Returns (nodes the gene sits on, where it starts in the first of them), or None when the gene
    does not sit on the path the way the context extraction will assume it does.
    """
    nodes_for_gene, start_in_first_node = geometry.nodes_covering(
        nodes_in_path, gene_contig_match.start - origin, gene_contig_match.end - origin)
    if not geometry.holds_gene(nodes_for_gene, start_in_first_node, gene_contig_match.aligned_length):
        return None
    return nodes_for_gene, start_in_first_node


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


def anchor_segment_in_contig(geometry, segment, placements_by_node, max_disagreement=MAX_ANCHOR_DISAGREEMENT):
    """Where a path segment of a gap-containing contig starts in contig coordinates.

    The gaps between segments are of unknown length, so a segment's origin cannot be derived from
    node lengths - it has to come from an aligned node. Any node of the segment will do: its offset
    within the segment is known exactly from the sequence walk, so the origin it implies is
    origin = where the node aligned - its offset in the segment. Using any node rather than only the
    first matters because short nodes never align (minimap2 reports nothing below ~200bp) and
    segments regularly start with one.

    Returns None when no node of the segment aligned, or when the candidate origins do not agree:
    a node that also aligned to a repeat elsewhere on the contig contributes a candidate far from
    the rest, and a segment whose candidates are spread out is not reliably anywhere. The median is
    taken rather than the best-scoring candidate because agreement, not alignment quality, is what
    says the segment is really there.
    """
    offsets = geometry.offsets(segment)
    origins = sorted(placement.origin - offset
                     for node, offset in zip(segment, offsets)
                     for placement in placements_by_node.get(node, []))
    if not origins:
        return None
    origin = origins[len(origins) // 2]
    agreeing = sum(1 for candidate in origins if abs(candidate - origin) <= max_disagreement)
    if agreeing * 2 <= len(origins):  # no strict majority agrees on where the segment starts
        log.debug(f'segment starting at {segment[0]} has alignments implying origins {origins}; '
                  f'too spread out to anchor it')
        return None
    return origin


def locate_gene_in_gappy_contig(geometry, segments, gene_contig_match, placements_by_node):
    """Locates a gene on a contig that SPAdes assembled from several graph paths joined using
    paired-end evidence.

    A segment is a run of nodes that really are adjacent in the graph - ';' in contigs.paths marks
    exactly the joins that are not graph edges - so requiring the gene's whole span to fall inside
    one segment is what makes its context a path the graph actually supports. A gene lying across a
    join belongs to no segment and is not located here; its context comes from the contig instead.
    """
    for segment in segments:
        origin = anchor_segment_in_contig(geometry, segment, placements_by_node)
        if origin is None:
            continue
        if origin <= gene_contig_match.start and gene_contig_match.end <= origin + geometry.length(segment):
            return locate_gene_on_path(geometry, segment, gene_contig_match, origin)
    return None


def add_node_list_to_genes_to_contigs(genes_to_contigs: Iterator[mc.GeneContigMatch],
                                      parsed_paths: Dict[str, List[List[str]]],
                                      geometry: pu.PathGeometry,
                                      placements_by_contig: Dict[str, Dict[str, list]],
                                      contigs_with_gaps: Set[str] = frozenset()):
    """Adds a nodes_list and a start_in_first_node to every gene-contig match that can be located in
    the assembly graph, and returns the matches worth analyzing.

    A match on one of contigs_with_gaps is kept even when it could not be located, because its
    context is taken from the contig rather than from the graph.
    """
    matches_to_analyze = []
    genes_located_per_route = Counter()

    # grouped by contig so that the geometry's memoized walks of a contig's path are reused by every
    # gene on it rather than recomputed per gene
    for contig, matches in groupby(sorted(genes_to_contigs, key=lambda match: match.contig),
                                   key=lambda match: match.contig):
        segments = parsed_paths.get(contig, [])
        contig_is_gappy = contig in contigs_with_gaps
        for gene_contig_match in matches:
            if contig_is_gappy:
                located = locate_gene_in_gappy_contig(geometry, segments, gene_contig_match,
                                                      placements_by_contig.get(contig, {}))
                route = 'a gap-free stretch of a gap-containing contig'
            elif segments:
                located = locate_gene_on_path(geometry, segments[0], gene_contig_match)
                route = 'the contig path'
            else:  # the contig has no entry in contigs.paths at all
                located = None
                route = None
            if located is not None:
                gene_contig_match.nodes_list, gene_contig_match.start_in_first_node = located
                genes_located_per_route[route] += 1
            if located is not None or contig_is_gappy:
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
    geometry = pu.PathGeometry(assembly_graph_nodes)
    # find nodes in contigs with gaps
    parsed_paths, contigs_with_gaps = pu.parse_paths_file(
        c.PATHS_PATH_TEMPLATE.format(assembly_dir=assembly_dir))
    placements_by_contig = map_nodes_to_contigs_w_gaps(contigs_with_gaps, assembly_graph_path, contigs_path, n_threads,
                                                       nodes_to_contigs_w_gaps_path)

    genes_to_analyze = add_node_list_to_genes_to_contigs(genes_to_contigs, parsed_paths, geometry,
                                                         placements_by_contig, contigs_with_gaps)
    n_located = sum(1 for match in genes_to_analyze if match.start_in_first_node is not None)
    log.info(f'found locations in the assembly graph for {n_located} genes, and kept '
             f'{len(genes_to_analyze) - n_located} more genes that were found on gap-containing contigs')
    return assembly_graph, genes_to_analyze, geometry, contigs_with_gaps


@pu.step_timing
def map_nodes_to_contigs_w_gaps(contigs_with_gaps, assembly_graph_path, contigs_path, n_threads,
                                nodes_to_contigs_w_gaps_path):
    """Aligns the graph's nodes to the gap-containing contigs, and returns where each node landed as
    {contig: {node: [NodePlacement]}}.

    Only gap-containing contigs are mapped, because they are the only ones whose path segments need
    anchoring - every other contig's path already places its nodes exactly.

    Nodes are named and oriented the way the rest of the pipeline names them, here rather than at
    every use: a placement's node is the one that runs in the contig's direction, and its origin is
    where the node itself starts, not where the alignment does.
    """
    # filter contigs fasta to keep only contigs with gaps
    contigs_w_gaps_path = contigs_path.replace('.fasta', '_w_gaps.fasta')
    SeqIO.write([contig for contig in SeqIO.parse(contigs_path, 'fasta') if contig.id in contigs_with_gaps],
                contigs_w_gaps_path, 'fasta')
    # find nodes in contigs with gaps
    nodes_to_contigs_path = sau.map_nodes_to_contigs(assembly_graph_path, contigs_w_gaps_path,
                                                     nodes_to_contigs_w_gaps_path,
                                                     nthreads=n_threads)
    return placements_from_minimap_results(pu.minimap_results_from_path(nodes_to_contigs_path))


def placements_from_minimap_results(nodes_to_contigs_df, score_threshold=NODE_TO_CONTIG_SCORE_TH):
    placements_by_contig = defaultdict(lambda: defaultdict(list))
    if not len(nodes_to_contigs_df):
        return placements_by_contig
    for record in nodes_to_contigs_df.itertuples():
        if record.mlen / record.qlen <= score_threshold:
            continue
        # every edge is in the fastg twice, as itself and as its reverse complement, and both align
        # here - they resolve to the same oriented node and become two placements of it in the same
        # spot, which anchor_segment_in_contig reads as two candidates that agree
        placement = mc.NodePlacement(node=node_oriented_with_contig(record.qname, record.strand),
                                     # minimap2 clips an alignment's ends, so the node itself starts
                                     # qstart before the alignment does
                                     origin=record.tstart - record.qstart,
                                     score=record.mlen / record.qlen)
        placements_by_contig[record.tname][placement.node].append(placement)
    return placements_by_contig
