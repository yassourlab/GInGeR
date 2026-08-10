import networkx as nx
from ginger import pipeline_utils as pu
import logging
from typing import Set
from Bio import SeqIO

log = logging.getLogger(__name__)

CONTIG_FALLBACK_PATH_NAME = 'contigfallback'
# the k-mer overlap two adjacent nodes are assumed to share while enumerating paths. The real overlap
# is found per pair by pipeline_utils.get_sequence_overlap, which is too expensive to call per step of
# the traversal
ASSUMED_NODE_OVERLAP = 55


def context_name(gene_contigs_match, path_name, side):
    """The name a context is written under - which is everything downstream knows about it, including
    the copy of the gene it was cut from (two copies can sit on the same nodes).

    '|'-separated with the gene first, so PathRefGenomeMatch can split it from the right whatever '|'
    the gene's own name contains. No other field may contain one, and nodes are joined with '_' rather
    than contigs.paths' ',' so a name needs no quoting in the output csvs.
    """
    match = gene_contigs_match
    return '|'.join([match.gene, match.contig, str(match.start), str(match.end), f'{match.score:.4f}',
                     '_'.join(match.nodes_list or []), path_name, side])


def save_paths_to_fasta_io_paths_approach(paths, paths_fasta_name, geometry, context_len,
                                          in_or_out, covered_by_gene, gene_contigs_match):
    """Writes the context_len bases of every path that flank the gene, and returns the names written.

    Contexts are always exactly context_len long, so a path with less than that left once the
    gene-covered part is trimmed off is dropped - fewer names come back than paths went in."""
    written_context_names = []

    with open(paths_fasta_name, 'a') as f:  # appended to: called once per gene location in the graph
        for path in paths:
            seq = geometry.sequence(path)
            covered_by_gene_int = int(covered_by_gene)
            # trim off the gene-covered portion, then take context_len bases off the end that flanks
            # the gene, so that the length check reflects the actual context that gets written
            if in_or_out == 'in':
                seq = seq[:len(seq) - covered_by_gene_int][-context_len:]
            else:
                seq = seq[covered_by_gene_int:][:context_len]
            if len(seq) == context_len:
                name = context_name(gene_contigs_match, '_'.join(path), in_or_out)
                written_context_names.append(name)
                f.write(f'>{name}\n{seq}\n')

    return written_context_names


def write_contexts_from_contig(contigs_index, gene_contigs_match, in_paths_fasta, out_paths_fasta, context_len):
    """Writes both of a gene's flanks, sliced straight out of the contig it was found on, and returns
    the 0, 1 or 2 names written - a side with less than context_len of flanking sequence is skipped.

    Done for genes on gap-containing contigs, where the graph describes the context poorly or not at
    all. Such a context may cross a paired-end-inferred join rather than a graph edge, which is weaker
    evidence - hence the identifiable 'contigfallback' path name.
    """
    contig_seq = str(contigs_index[gene_contigs_match.contig].seq)
    sides = [('in', in_paths_fasta,
              contig_seq[max(0, gene_contigs_match.start - context_len):gene_contigs_match.start]),
             ('out', out_paths_fasta,
              contig_seq[gene_contigs_match.end:gene_contigs_match.end + context_len])]

    written = []
    for side, paths_fasta, seq in sides:
        if len(seq) == context_len:
            name = context_name(gene_contigs_match, CONTIG_FALLBACK_PATH_NAME, side)
            with open(paths_fasta, 'a') as f:
                f.write(f'>{name}\n{seq}\n')
            written.append(name)
    return written


def paths_enumerator(graph, stack, max_depth, max_length, neighbors_func, reverse=False, covered_by_gene=0):
    out_paths = []
    while stack:
        (vertex, path) = stack.pop()
        # an estimate, not the length: it assumes every join collapses ASSUMED_NODE_OVERLAP bases and
        # charges one join per node rather than per pair, so it runs a little short. Only the bound on
        # how far to keep walking depends on it - what gets written is cut to length from the real
        # sequence in save_paths_to_fasta_io_paths_approach
        total_length_estimation = sum([length for node, length in path]) - (
                len(path) * ASSUMED_NODE_OVERLAP) - covered_by_gene
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
def extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, geometry: pu.PathGeometry,
                                                      genes_to_contigs, depth_limit, context_len,
                                                      in_paths_fasta, out_paths_fasta, contigs_path,
                                                      contigs_with_gaps: Set[str] = frozenset()):
    """Writes a context candidate fasta for the incoming and for the outgoing side of every gene, and
    returns the length of every gene.

    Every context is exactly context_len long, and a side that cannot supply that much gets no context.
    Contexts are read off the assembly graph; a gene on one of contigs_with_gaps additionally gets
    contexts sliced out of the contig (see write_contexts_from_contig) - pass an empty set to turn that
    off.

    Every copy of a gene is handled on its own even when two copies sit on the same nodes: pairing one
    copy's incoming context with another's outgoing context would describe a stretch that is in no
    contig.
    """
    gene_lengths = {}
    n_contig_contexts = n_graph_contexts = 0
    contigs_index = SeqIO.index(contigs_path, 'fasta') if contigs_with_gaps else None
    # start from empty fastas - everything below appends to them
    for fasta_file in [in_paths_fasta, out_paths_fasta]:
        open(fasta_file, 'w').close()
    try:
        for gene_contigs_match in genes_to_contigs:
            gene_lengths[gene_contigs_match.gene] = gene_contigs_match.gene_length
            if gene_contigs_match.contig in contigs_with_gaps:
                n_contig_contexts += len(write_contexts_from_contig(contigs_index, gene_contigs_match,
                                                                    in_paths_fasta, out_paths_fasta, context_len))

            if gene_contigs_match.start_in_first_node is None:  # the gene was not located in the graph
                continue

            # unpacking variables
            first_node = gene_contigs_match.nodes_list[0]
            last_node = gene_contigs_match.nodes_list[-1]
            start_in_first_node = gene_contigs_match.start_in_first_node
            gene_nodes_length = geometry.length(gene_contigs_match.nodes_list)

            # in paths
            in_covered_by_gene = assembly_graph.nodes[first_node]['length'] - start_in_first_node
            in_paths_initial_stack = [(first_node, [(first_node, assembly_graph.nodes[first_node]['length'])])]
            in_paths = paths_enumerator(assembly_graph, in_paths_initial_stack, depth_limit, context_len,
                                        nx.DiGraph.predecessors, reverse=True, covered_by_gene=in_covered_by_gene)
            n_graph_contexts += len(save_paths_to_fasta_io_paths_approach(
                in_paths, in_paths_fasta, geometry, context_len, 'in', in_covered_by_gene, gene_contigs_match))

            # out paths
            out_paths_initial_stack = [(last_node, [(last_node, assembly_graph.nodes[last_node]['length'])])]
            # the outgoing context has to start where the gene's alignment on the contig ends, which
            # is aligned_length - not gene_length - away from its start. gene_length is the length of
            # the reference protein and can be a good deal longer than the aligned span, which would
            # start the context that many bases past the end of the gene and drop them from the
            # sequence altogether.
            # PathGeometry.holds_gene checked when the gene was located that it ends on these nodes
            # and reaches the last of them, so this is between 0 and that node's length
            bases_after_gene_in_path = gene_nodes_length - (start_in_first_node +
                                                            gene_contigs_match.aligned_length)
            out_covered_by_gene = assembly_graph.nodes[last_node]['length'] - bases_after_gene_in_path
            out_paths = paths_enumerator(assembly_graph, out_paths_initial_stack, depth_limit, context_len,
                                         nx.DiGraph.successors, covered_by_gene=out_covered_by_gene)
            n_graph_contexts += len(save_paths_to_fasta_io_paths_approach(
                out_paths, out_paths_fasta, geometry, context_len, 'out', out_covered_by_gene, gene_contigs_match))

            if len(in_paths) == 0 or len(out_paths) == 0:
                log.info(f'{gene_contigs_match} in {len(in_paths)} out {len(out_paths)}')
    finally:
        if contigs_index is not None:
            contigs_index.close()

    log.info(f'read {n_graph_contexts} contexts off the assembly graph and took {n_contig_contexts} more from the '
             f'sequence of gap-containing contigs')
    return gene_lengths
