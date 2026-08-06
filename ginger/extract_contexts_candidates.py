import networkx as nx
from ginger import pipeline_utils as pu
import logging
import datetime as dt
from typing import Dict, Set
from Bio import SeqIO

log = logging.getLogger(__name__)


def save_paths_to_fasta_io_paths_approach(paths, paths_fasta_name, records_dict, context_len, node_to_find=None,
                                          in_or_out=None, covered_by_gene=0, gene_and_node='', match_score=None):
    """Writes the context_len bases of every path that flank the gene. Contexts are always exactly
    context_len long - a path that has less than that left once the gene-covered part is trimmed off
    is dropped."""
    in_paths_lengths = {}
    node_locations = {}

    with open(paths_fasta_name, 'a') as f:  # there is an 'a' here because I call this function once per gene location in the graph
        for path in paths:
            seq, node_start = pu.generate_str_from_list_of_nodes(records_dict, path, node_to_find)
            covered_by_gene_int = int(covered_by_gene)
            # trim off the gene-covered portion, then take context_len bases off the end that flanks
            # the gene, so that the length check reflects the actual context that gets written
            if in_or_out == 'in':
                seq = seq[:len(seq) - covered_by_gene_int][-context_len:]
            if in_or_out == 'out':
                seq = seq[covered_by_gene_int:][:context_len]
            if len(seq) == context_len:
                node_locations['_'.join(path)] = node_start
                in_paths_lengths['_'.join(path)] = len(seq)
                match_score_str = f"_match_{match_score:.4f}" if match_score is not None else ""
                f.write(f">{gene_and_node}{match_score_str}_path_{'_'.join(path)}\n" if gene_and_node else f">{'_'.join(path)}\n")
                f.write(f'{seq}\n')

    return in_paths_lengths, node_locations


def gene_and_nodes_name(gene_contigs_match):
    """The '{gene}_nodes_{nodes}' prefix of a context's name in the output fasta. A gene that could
    not be located in the assembly graph has no nodes, and gets an empty nodes part."""
    return f"{gene_contigs_match.gene}_nodes_{'_'.join(gene_contigs_match.nodes_list or [])}"


def write_contexts_from_contig(contigs_index, gene_contigs_match, in_paths_fasta, out_paths_fasta, context_len):
    """Writes both of a gene's flanks, sliced straight out of the contig it was found on.

    This is done for every gene found on a gap-containing contig - a contig SPAdes assembled from
    several graph paths joined using paired-end evidence. On such a contig the gene sits on a node
    that is usually a dead end in the graph, and sometimes it can't be located in the graph at all,
    so the graph describes its context poorly or not at all, while the contig has flanking sequence
    on both sides. A context taken from the contig may cross one of those paired-end-inferred joins,
    which is weaker evidence than pure graph sequence - hence the identifiable 'contigfallback' path
    name it gets in the output.

    Returns how many contexts were written (0, 1 or 2 - a side that has less than context_len of
    flanking sequence left in the contig is skipped).
    """
    contig_seq = str(contigs_index[gene_contigs_match.contig].seq)
    name = (f'{gene_and_nodes_name(gene_contigs_match)}_match_{gene_contigs_match.score:.4f}_path_'
            f'contigfallback_{gene_contigs_match.contig}_{gene_contigs_match.start}_{gene_contigs_match.end}')
    sides = [(in_paths_fasta, contig_seq[max(0, gene_contigs_match.start - context_len):gene_contigs_match.start]),
             (out_paths_fasta, contig_seq[gene_contigs_match.end:gene_contigs_match.end + context_len])]

    n_written = 0
    for paths_fasta, seq in sides:
        if len(seq) == context_len:
            with open(paths_fasta, 'a') as f:
                f.write(f'>{name}\n{seq}\n')
            n_written += 1
    return n_written


def paths_enumerator(graph, stack, max_depth, max_length, neighbors_func, reverse=False, covered_by_gene=0):
    out_paths = []
    while stack:
        (vertex, path) = stack.pop()
        total_length_estimation = sum([length for node, length in path]) - (len(path) * 55) - covered_by_gene
        neighbors = list(neighbors_func(graph, vertex))
        if len(neighbors) > 0 and total_length_estimation < max_length and len(path) < max_depth:
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
                                                      genes_to_contigs, depth_limit, context_len,
                                                      in_paths_fasta, out_paths_fasta, contigs_path,
                                                      contigs_with_gaps: Set[str] = frozenset()):
    """Writes a context candidate fasta for the incoming and for the outgoing side of every gene.

    Every context is exactly context_len long, and a side that cannot supply that much sequence gets
    no context. Contexts are read off the assembly graph. A gene found on one of contigs_with_gaps additionally
    gets contexts sliced straight out of the contig, whether or not the graph could locate it or
    describe its context (see write_contexts_from_contig). Pass an empty contigs_with_gaps to turn
    that off.
    """
    gene_and_nodes_path_set = set()
    gene_lengths = {}
    n_contig_contexts = 0
    contigs_index = SeqIO.index(contigs_path, 'fasta') if contigs_with_gaps else None
    # start from empty fastas - everything below appends to them
    for fasta_file in [in_paths_fasta, out_paths_fasta]:
        open(fasta_file, 'w').close()
    try:
        for gene_contigs_match in genes_to_contigs:
            gene_lengths[gene_contigs_match.gene] = gene_contigs_match.gene_length
            if gene_contigs_match.contig in contigs_with_gaps:
                n_contig_contexts += write_contexts_from_contig(contigs_index, gene_contigs_match, in_paths_fasta,
                                                                out_paths_fasta, context_len)

            gene_and_nodes_path_str = gene_and_nodes_name(gene_contigs_match)
            if gene_contigs_match.start_in_first_node is None:  # the gene was not located in the graph
                log.info(f'{dt.datetime.now()} did not extract contexts from the graph for '
                         f'{gene_and_nodes_path_str} because start_in_first_node is None')
                continue
            if gene_and_nodes_path_str in gene_and_nodes_path_set:  # already ran the pipeline for this gene location
                continue
            gene_and_nodes_path_set.add(gene_and_nodes_path_str)

            # unpacking variables
            first_node = gene_contigs_match.nodes_list[0]
            last_node = gene_contigs_match.nodes_list[-1]
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
            save_paths_to_fasta_io_paths_approach(in_paths, in_paths_fasta, nodes_with_edges_and_sequences, context_len,
                                                  in_or_out='in', covered_by_gene=in_covered_by_gene,
                                                  gene_and_node=gene_and_nodes_path_str,
                                                  match_score=gene_contigs_match.score)

            # out paths
            out_paths_initial_stack = [(last_node, [(last_node, assembly_graph.nodes[last_node]['length'])])]
            gene_end_in_last_node = gene_nodes_length - start_in_first_node - gene_length
            out_covered_by_gene = assembly_graph.nodes[last_node]['length'] - gene_end_in_last_node
            out_paths = paths_enumerator(assembly_graph, out_paths_initial_stack, depth_limit, context_len,
                                         nx.DiGraph.successors, covered_by_gene=out_covered_by_gene)
            save_paths_to_fasta_io_paths_approach(out_paths, out_paths_fasta, nodes_with_edges_and_sequences, context_len,
                                                  in_or_out='out', covered_by_gene=out_covered_by_gene,
                                                  gene_and_node=gene_and_nodes_path_str,
                                                  match_score=gene_contigs_match.score)

            if len(in_paths) == 0 or len(out_paths) == 0:
                log.info(f'{gene_contigs_match} in {len(in_paths)} out {len(out_paths)}')
    finally:
        if contigs_index is not None:
            contigs_index.close()

    log.info(f'took {n_contig_contexts} contexts from the sequence of gap-containing contigs')
    return gene_lengths
