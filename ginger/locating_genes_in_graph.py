import pyfastg
from ginger import pipeline_utils as pu
import logging
from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import constants as c
from Bio import SeqIO
from typing import Dict, Iterator, List
import networkx as nx
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


def add_location_in_graph_based_on_contigs_paths(nodes_sequences_dict, nodes_in_path, gene_contig_match):
    #  this is not the exact start but it's good enough
    contig_start = gene_contig_match.start
    contig_end = gene_contig_match.end
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


def format_contig_name_for_nodes_list(contig_name):
    contig_num = contig_name.split('_')[1]
    if contig_num[-1] == "'":
        return contig_num + '-'
    else:
        return contig_num + '+'


def get_start_in_first_node(gene_contig_match, top_nodes_to_contigs_match):
    if top_nodes_to_contigs_match.strand == '+':
        return gene_contig_match.start - top_nodes_to_contigs_match.contig_start
    else:
        return top_nodes_to_contigs_match.contig_end - gene_contig_match.end


def add_location_in_graph_based_on_nodes_to_contigs(gene_contig_match, nodes_to_contigs_df):
    pos_stran_conditions = ((nodes_to_contigs_df.contig_start <= gene_contig_match.start) & (
            gene_contig_match.start <= nodes_to_contigs_df.contig_end) & (nodes_to_contigs_df.strand == '+'))
    neg_strand_conditions = ((nodes_to_contigs_df.contig_start <= gene_contig_match.end) & (
            gene_contig_match.end <= nodes_to_contigs_df.contig_end) & (nodes_to_contigs_df.strand == '-'))
    or_between_conditions = pos_stran_conditions | neg_strand_conditions
    filtered_nodes_to_contigs = nodes_to_contigs_df[
        (nodes_to_contigs_df.contig == gene_contig_match.contig) & or_between_conditions]

    if filtered_nodes_to_contigs.empty:
        return None
    else:
        top_nodes_to_contigs_match = filtered_nodes_to_contigs.iloc[0]
        gene_contig_match.nodes_list = [format_contig_name_for_nodes_list(top_nodes_to_contigs_match.node)]
        gene_contig_match.start_in_first_node = get_start_in_first_node(gene_contig_match, top_nodes_to_contigs_match)
    return gene_contig_match


def add_node_list_to_genes_to_contigs(genes_to_contigs: Iterator[mc.GeneContigMatch],
                                      parsed_paths_without_gaps: Dict[str, List],
                                      nodes_sequences_dict: Dict[str, SeqIO.SeqRecord],
                                      nodes_to_contigs_df: pd.DataFrame):
    matches_with_nodes_list_and_start_location = []

    for gene_contig_match in genes_to_contigs:
        contig_path_in_graph = parsed_paths_without_gaps.get(gene_contig_match.contig, None)
        if contig_path_in_graph:
            updated_gene_contig_match = add_location_in_graph_based_on_contigs_paths(nodes_sequences_dict,
                                                                                     contig_path_in_graph,
                                                                                     gene_contig_match)
        else:
            updated_gene_contig_match = add_location_in_graph_based_on_nodes_to_contigs(gene_contig_match,
                                                                                        nodes_to_contigs_df)
        if updated_gene_contig_match is not None:
            matches_with_nodes_list_and_start_location.append(updated_gene_contig_match)
        # else:
        #     log.warning(f'gene {gene_contig_match.gene} not found in the assembly graph')
    return matches_with_nodes_list_and_start_location


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
                          pident_filtering_th: float, genes_to_contigs_path: str) -> Iterator[mc.GeneContigMatch]:
    if not os.path.exists(genes_to_contigs_path) or not os.path.isfile(genes_to_contigs_path):
        log.info(f'running mmseqs2 to find genes in contigs')
        genes_to_contigs_path = sau.map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path,
                                                         nthreads=n_minimap_threads)
    genes_to_contigs = sau.read_and_filter_mmseq2_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                          pident_filtering_th)
    return genes_to_contigs


@pu.step_timing
def locate_genes_in_graph(assembly_dir: str, gene_pident_filtering_th: float, genes_path: str, n_threads: int,
                          temp_folder: str):  # -> Tuple[networkx.DiGraph,??? ,Dict[str, SeqIO.SeqRecord]]
    contigs_path = c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    assembly_graph_path = c.ASSEMBLY_GRAPH_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    genes_to_contigs_path = c.GENES_TO_CONTIGS_TEMPLATE.format(temp_files_path=temp_folder)
    nodes_to_contigs_w_gaps_path = c.NODES_TO_CONTIGS_W_GAPS_TEMPLATE.format(temp_files_path=temp_folder)

    # find genes in contigs
    genes_to_contigs = find_genes_in_contigs(genes_path, contigs_path, n_threads,
                                             gene_pident_filtering_th, genes_to_contigs_path)

    log.info(f'found {len(set([m.gene for m in genes_to_contigs]))} genes in the assembly graph')
    if genes_to_contigs is None:
        return None, None, None

    assembly_graph_nodes = get_nodes_dict_from_fastg_file(assembly_graph_path)
    assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    # find nodes in contigs with gaps
    parsed_paths_without_gaps, contigs_with_gaps = pu.parse_paths_file(
        c.PATHS_PATH_TEMPLATE.format(assembly_dir=assembly_dir), assembly_graph.nodes)
    nodes_to_contigs_df = map_nodes_to_contigs_w_gaps(contigs_with_gaps, assembly_graph_path, contigs_path, n_threads,
                                                      nodes_to_contigs_w_gaps_path)

    genes_with_location_in_graph = add_node_list_to_genes_to_contigs(genes_to_contigs, parsed_paths_without_gaps,
                                                                     assembly_graph_nodes, nodes_to_contigs_df)
    log.info(f'found locations in the assembly graph for {len(genes_with_location_in_graph)} genes')
    return assembly_graph, genes_with_location_in_graph, assembly_graph_nodes


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
    nodes_to_contigs_df['node'] = nodes_to_contigs_df.qname.apply(lambda x: x[:-1].split(':')[0])
    nodes_to_contigs_df['score'] = nodes_to_contigs_df.mlen / nodes_to_contigs_df.qlen
    nodes_to_contigs_df = \
        nodes_to_contigs_df[(nodes_to_contigs_df.score > 0.95)].rename(
            columns={'tname': 'contig', 'tstart': 'contig_start', 'tend': 'contig_end'})[
            ['contig', 'contig_start', 'contig_end', 'node', 'score', 'strand']]
    nodes_to_contigs_df.sort_values(['score', 'node'], inplace=True, ascending=(False, True))
    nodes_to_contigs_df.drop_duplicates(subset=['contig', 'contig_start', 'contig_end', 'score'],
                                        keep='first', inplace=True)
    return nodes_to_contigs_df
