import pyfastg
import pipeline_utils as pu
import logging
import sequence_alignment_utils as sau
import matches_classes as mc
import constants as c
from Bio import SeqIO
from typing import Dict, Iterator, Tuple
import networkx as nx

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


def add_nodes_list_and_start_location_to_gene_contig_match(nodes_sequences_dict, nodes_in_path, gene_contig_match):
    #  this is not the exact start but it's good enough
    contig_start = gene_contig_match.contig_start
    contig_end = gene_contig_match.contig_end
    if not nodes_in_path:
        return gene_contig_match
    elif len(nodes_in_path) == 1:
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
    # TODO I'm modifying and then returning the same object. I think it's not optimal
    return gene_contig_match


def get_genes_to_contigs_with_nodes_list(genes_to_contigs: Iterator[mc.GeneContigMatch], assembly_graph: nx.DiGraph,
                                         nodes_sequences_dict: Dict[str, SeqIO.SeqRecord], assembly_dir: str):
    contig_nodes_dict = pu.get_contig_nodes_dict(assembly_graph.nodes,
                                                 c.PATHS_PATH_TEMPLATE.format(assembly_dir=assembly_dir),
                                                 keep_contigs_with_gaps=False)
    matches_with_nodes_list_and_start_location = []
    for gene_contig_match in genes_to_contigs:  # TODO turn this to an iterator when I finish debugging
        updated_gene_contig_match = add_nodes_list_and_start_location_to_gene_contig_match(nodes_sequences_dict,
                                                                                           contig_nodes_dict.get(
                                                                                               gene_contig_match.contig,
                                                                                               None), gene_contig_match)
        if updated_gene_contig_match.nodes_list and updated_gene_contig_match.start_in_first_node:
            matches_with_nodes_list_and_start_location.append(updated_gene_contig_match)
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


def find_genes_in_contigs(temp_dir: str, genes_path: str, contigs_path: str, n_minimap_threads: int,
                          pident_filtering_th: float) -> Iterator[mc.GeneContigMatch]:
    genes_to_contigs_path = c.GENES_TO_CONTIGS_TEMPLATE.format(temp_files_path=temp_dir)

    genes_to_contigs_path = sau.map_genes_to_contexts(genes_path, contigs_path, genes_to_contigs_path,
                                                      nthreads=n_minimap_threads)
    genes_to_contigs = sau.read_and_filter_minimap_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                           pident_filtering_th,
                                                           verbose=True)
    return genes_to_contigs


def locate_genes_in_graph(assembly_dir: str, gene_pident_filtering_th: float, genes_path: str, n_minimap_threads: int,
                          temp_folder: str):  # -> Tuple[networkx.DiGraph,??? ,???]
    contigs_path = c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    assembly_graph_path = c.ASSEMBLY_GRAPH_PATH_TEMPLATE.format(assembly_dir=assembly_dir)
    nodes_with_edges_and_sequences = get_nodes_dict_from_fastg_file(assembly_graph_path)

    genes_to_contigs = find_genes_in_contigs(temp_folder, genes_path, contigs_path, n_minimap_threads,
                                             gene_pident_filtering_th)

    assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    genes_with_location_in_graph = get_genes_to_contigs_with_nodes_list(genes_to_contigs, assembly_graph,
                                                                        nodes_with_edges_and_sequences, assembly_dir)
    return assembly_graph, genes_with_location_in_graph, nodes_with_edges_and_sequences
