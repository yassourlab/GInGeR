import networkx as nx
from ginger import pipeline_utils as pu
from ginger import matches_classes as mc
import logging
import datetime as dt
from collections import namedtuple
from typing import Dict, Set
from Bio import SeqIO

log = logging.getLogger(__name__)

CONTEXTS_TO_LOCI_COLUMNS = ['context_name', 'side', 'contig', 'gene_start', 'gene_end', 'nodes_list',
                            'match_score', 'source']
ContextLocusRow = namedtuple('ContextLocusRow', CONTEXTS_TO_LOCI_COLUMNS)


def context_locus_rows(written_contexts, gene_contigs_match, source):
    """One row per written context, recording the gene copy it was cut from.

    A context's name says which gene and which graph nodes it came from, but that is not enough to
    identify the copy: a gene that could not be located in the graph has no nodes in its name at
    all, and two copies can sit on the same nodes. Only the contig interval identifies it, and only
    here is it still known - so it is written down rather than parsed back out of the name later.
    """
    return [ContextLocusRow(context_name, side, gene_contigs_match.contig, gene_contigs_match.start,
                            gene_contigs_match.end, '_'.join(gene_contigs_match.nodes_list or []),
                            gene_contigs_match.score, source)
            for side, context_name in written_contexts]


def write_contexts_to_loci_table(rows, contexts_to_loci_path):
    with open(contexts_to_loci_path, 'w') as f:
        f.write('\t'.join(CONTEXTS_TO_LOCI_COLUMNS) + '\n')
        for row in rows:
            f.write('\t'.join(str(field) for field in row) + '\n')


def save_paths_to_fasta_io_paths_approach(paths, paths_fasta_name, records_dict, context_len, node_to_find=None,
                                          in_or_out=None, covered_by_gene=0, gene_and_node='', match_score=None):
    """Writes the context_len bases of every path that flank the gene. Contexts are always exactly
    context_len long - a path that has less than that left once the gene-covered part is trimmed off
    is dropped.

    The lengths are keyed by the name the context was written under, so that a caller can record
    what it wrote without having to rebuild those names itself."""
    written_context_lengths = {}
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
                path_name = '_'.join(path)
                match_score_str = f"_match_{match_score:.4f}" if match_score is not None else ""
                context_name = f'{gene_and_node}{match_score_str}_path_{path_name}' if gene_and_node else path_name
                node_locations[path_name] = node_start
                written_context_lengths[context_name] = len(seq)
                f.write(f'>{context_name}\n{seq}\n')

    return written_context_lengths, node_locations


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

    Returns the (side, context name) pairs written - 0, 1 or 2 of them, since a side that has less
    than context_len of flanking sequence left in the contig is skipped. Both sides are written
    under the same name, as they come from the same place.
    """
    contig_seq = str(contigs_index[gene_contigs_match.contig].seq)
    name = (f'{gene_and_nodes_name(gene_contigs_match)}_match_{gene_contigs_match.score:.4f}_path_'
            f'contigfallback_{gene_contigs_match.contig}_{gene_contigs_match.start}_{gene_contigs_match.end}')
    sides = [('in', in_paths_fasta,
              contig_seq[max(0, gene_contigs_match.start - context_len):gene_contigs_match.start]),
             ('out', out_paths_fasta,
              contig_seq[gene_contigs_match.end:gene_contigs_match.end + context_len])]

    written = []
    for side, paths_fasta, seq in sides:
        if len(seq) == context_len:
            with open(paths_fasta, 'a') as f:
                f.write(f'>{name}\n{seq}\n')
            written.append((side, name))
    return written


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
                                                      contigs_with_gaps: Set[str] = frozenset(),
                                                      contexts_to_loci_path=None):
    """Writes a context candidate fasta for the incoming and for the outgoing side of every gene.

    Every context is exactly context_len long, and a side that cannot supply that much sequence gets
    no context. Contexts are read off the assembly graph. A gene found on one of contigs_with_gaps additionally
    gets contexts sliced straight out of the contig, whether or not the graph could locate it or
    describe its context (see write_contexts_from_contig). Pass an empty contigs_with_gaps to turn
    that off.

    Returns the length of every gene, and the locus every context was cut from - written to
    contexts_to_loci_path as well, when one is given.
    """
    gene_and_nodes_path_set = set()
    gene_lengths = {}
    locus_rows = []
    n_contig_contexts = 0
    contigs_index = SeqIO.index(contigs_path, 'fasta') if contigs_with_gaps else None
    # start from empty fastas - everything below appends to them
    for fasta_file in [in_paths_fasta, out_paths_fasta]:
        open(fasta_file, 'w').close()
    try:
        for gene_contigs_match in genes_to_contigs:
            gene_lengths[gene_contigs_match.gene] = gene_contigs_match.gene_length
            if gene_contigs_match.contig in contigs_with_gaps:
                written = write_contexts_from_contig(contigs_index, gene_contigs_match, in_paths_fasta,
                                                     out_paths_fasta, context_len)
                locus_rows.extend(context_locus_rows(written, gene_contigs_match, 'contigfallback'))
                n_contig_contexts += len(written)

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
            gene_nodes_length = len(
                pu.generate_str_from_list_of_nodes(nodes_with_edges_and_sequences, nodes_list_for_gene, None)[0])

            # in paths
            in_covered_by_gene = assembly_graph.nodes[first_node]['length'] - start_in_first_node
            in_paths_initial_stack = [(first_node, [(first_node, assembly_graph.nodes[first_node]['length'])])]
            in_paths = paths_enumerator(assembly_graph, in_paths_initial_stack, depth_limit, context_len,
                                        nx.DiGraph.predecessors, reverse=True, covered_by_gene=in_covered_by_gene)
            written_in_contexts, _ = save_paths_to_fasta_io_paths_approach(
                in_paths, in_paths_fasta, nodes_with_edges_and_sequences, context_len,
                in_or_out='in', covered_by_gene=in_covered_by_gene,
                gene_and_node=gene_and_nodes_path_str, match_score=gene_contigs_match.score)
            locus_rows.extend(context_locus_rows([('in', name) for name in written_in_contexts],
                                                 gene_contigs_match, 'graph'))

            # out paths
            out_paths_initial_stack = [(last_node, [(last_node, assembly_graph.nodes[last_node]['length'])])]
            # the outgoing context has to start where the gene's alignment on the contig ends, which
            # is aligned_length - not gene_length - away from its start. gene_length is the length of
            # the reference protein and can be a good deal longer than the aligned span, which would
            # start the context that many bases past the end of the gene and drop them from the
            # sequence altogether.
            bases_after_gene_in_path = gene_nodes_length - (start_in_first_node +
                                                            gene_contigs_match.aligned_length)
            if bases_after_gene_in_path < 0:
                log.info(f'{dt.datetime.now()} did not extract an outgoing context from the graph for '
                         f'{gene_and_nodes_path_str} because the gene alignment ends past the nodes it '
                         f'was located on')
                continue
            out_covered_by_gene = assembly_graph.nodes[last_node]['length'] - bases_after_gene_in_path
            out_paths = paths_enumerator(assembly_graph, out_paths_initial_stack, depth_limit, context_len,
                                         nx.DiGraph.successors, covered_by_gene=out_covered_by_gene)
            written_out_contexts, _ = save_paths_to_fasta_io_paths_approach(
                out_paths, out_paths_fasta, nodes_with_edges_and_sequences, context_len,
                in_or_out='out', covered_by_gene=out_covered_by_gene,
                gene_and_node=gene_and_nodes_path_str, match_score=gene_contigs_match.score)
            locus_rows.extend(context_locus_rows([('out', name) for name in written_out_contexts],
                                                 gene_contigs_match, 'graph'))

            if len(in_paths) == 0 or len(out_paths) == 0:
                log.info(f'{gene_contigs_match} in {len(in_paths)} out {len(out_paths)}')
    finally:
        if contigs_index is not None:
            contigs_index.close()

    log.info(f'took {n_contig_contexts} contexts from the sequence of gap-containing contigs')
    if contexts_to_loci_path is not None:
        write_contexts_to_loci_table(locus_rows, contexts_to_loci_path)
    contexts_to_loci = {row.context_name: mc.GeneLocus(row.contig, row.gene_start, row.gene_end)
                        for row in locus_rows}
    return gene_lengths, contexts_to_loci
