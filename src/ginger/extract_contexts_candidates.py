import networkx as nx
import pipeline_utils as pu
import math
import logging
import datetime as dt
from typing import Dict
from Bio import SeqIO

log = logging.getLogger(__name__)


def save_paths_to_fasta_io_paths_approach(paths, paths_fasta_name, records_dict, node_to_find=None,
                                          context_len=math.inf, in_or_out=None, covered_by_gene=0, gene_and_node=''):
    # TODO I think I don't really need the outputs of this function anymore
    # print(f'{dt.datetime.now()} saving paths to {paths_fasta_name}')
    in_paths_lengths = {}
    node_locations = {}
    with open(paths_fasta_name, 'a') as f:
        for path in paths:
            seq, node_start = pu.generate_str_from_list_of_nodes(records_dict, path, node_to_find)
            node_locations['_'.join(path)] = node_start
            if len(seq) > context_len:
                if in_or_out == 'in':
                    seq = seq[-(int(context_len + covered_by_gene)):-int(covered_by_gene)]
                if in_or_out == 'out':
                    seq = seq[int(covered_by_gene):int(covered_by_gene + context_len)]
            else:
                if in_or_out == 'in':
                    seq = seq[:-int(covered_by_gene)]
                if in_or_out == 'out':
                    seq = seq[int(covered_by_gene):]
            in_paths_lengths['_'.join(path)] = len(seq)
            f.write(f">{gene_and_node}_path_{'_'.join(path)}\n" if gene_and_node else f">{'_'.join(path)}\n")
            f.write(f'{seq}\n')

    return in_paths_lengths, node_locations


def paths_enumerator(graph, stack, max_depth, min_length, neighbors_func, reverse=False, covered_by_gene=0):
    # TODO in cases where there are no incoming or out coming paths, I think that we don't take the path that consists only of the node itself into account - for example gene_Subject_11486372__Scaffold_30890__Start_212__End_556_nodes_722037+ in the graph with just 1M short reads
    out_paths = []
    while stack:
        (vertex, path) = stack.pop()
        total_length_estimation = sum([length for node, length in path]) - (len(path) * 55) - covered_by_gene
        neighbors = list(neighbors_func(graph, vertex))
        if len(neighbors) > 0 and total_length_estimation < min_length and len(path) < max_depth:
            for neighbor in neighbors:
                if reverse:
                    stack.append((neighbor, [(neighbor, graph.nodes[neighbor]['length'])] + path))
                else:
                    stack.append((neighbor, path + [(neighbor, graph.nodes[neighbor]['length'])]))
        else:
            out_paths.append([n for n, l in path])
    return out_paths


@pu.step_timing
def extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                      nodes_with_edges_and_sequences: Dict[str, SeqIO.SeqRecord],
                                                      genes_to_contigs, depth_limit,
                                                      context_len,
                                                      in_paths_fasta, out_paths_fasta):
    gene_and_nodes_path_set = set()
    gene_lengths = {}
    for gene_contigs_match in genes_to_contigs:
        gene_and_nodes_path_str = f"{gene_contigs_match.gene}_nodes_{'_'.join(gene_contigs_match.nodes_list)}"
        if gene_and_nodes_path_str in gene_and_nodes_path_set or gene_contigs_match.start_in_first_node is None:  # I already ran the pipeline for this and there is no need to do it again
            log.debug(
                f'{dt.datetime.now()} pipeline did not run for {gene_and_nodes_path_str} {gene_contigs_match.start_in_first_node}')
        else:
            # unpacking variables
            gene_and_nodes_path_set.add(gene_and_nodes_path_str)
            first_node = gene_contigs_match.nodes_list[0]
            last_node = gene_contigs_match.nodes_list[-1]
            gene_name = gene_contigs_match.gene
            nodes_list_for_gene = gene_contigs_match.nodes_list
            start_in_first_node = gene_contigs_match.start_in_first_node
            gene_length = gene_contigs_match.gene_length
            gene_nodes_length = len(
                pu.generate_str_from_list_of_nodes(nodes_with_edges_and_sequences, nodes_list_for_gene, None)[0])

            # in paths
            in_covered_by_gene = assembly_graph.nodes[first_node]['length'] - start_in_first_node
            in_paths_initial_stack = [(first_node, [(first_node, assembly_graph.nodes[first_node]['length'])])]
            in_paths = paths_enumerator(assembly_graph, in_paths_initial_stack, depth_limit, context_len,
                                        nx.DiGraph.predecessors, reverse=True, covered_by_gene=in_covered_by_gene)
            in_paths_lengths, _ = save_paths_to_fasta_io_paths_approach(in_paths, in_paths_fasta,
                                                                        nodes_with_edges_and_sequences,
                                                                        context_len=context_len,
                                                                        in_or_out='in',
                                                                        covered_by_gene=in_covered_by_gene,
                                                                        gene_and_node=gene_and_nodes_path_str)

            # out paths
            out_paths_initial_stack = [(last_node, [(last_node, assembly_graph.nodes[last_node][
                'length'])])]
            gene_end_in_last_node = gene_nodes_length - start_in_first_node - gene_length
            out_covered_by_gene = assembly_graph.nodes[last_node]['length'] - gene_end_in_last_node

            out_paths = paths_enumerator(assembly_graph, out_paths_initial_stack, depth_limit, context_len,
                                         nx.DiGraph.successors, covered_by_gene=out_covered_by_gene)
            out_paths_lengths, _ = save_paths_to_fasta_io_paths_approach(out_paths, out_paths_fasta,
                                                                         nodes_with_edges_and_sequences,
                                                                         context_len=context_len, in_or_out='out',
                                                                         covered_by_gene=out_covered_by_gene,
                                                                         gene_and_node=gene_and_nodes_path_str)
            gene_lengths[gene_name] = gene_length
            if len(in_paths) == 0 or len(out_paths) == 0:
                log.debug(f'{gene_contigs_match} in {len(in_paths)} out {len(out_paths)}')

    return gene_lengths
