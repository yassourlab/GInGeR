import pyfastg
import pipeline_utils as pu
import datetime as dt
import logging
import sequence_alignment_utils as sau
import matches_classes as mc

log = logging.getLogger(__name__)


def add_nodes_list_and_start_location_to_gene_contig_match(records_dict, nodes_in_path, gene_contig_match):
    #  this is not the exact start but it's good enough
    contig_start = gene_contig_match.contig_start
    contig_end = gene_contig_match.contig_end
    if not nodes_in_path:
        return gene_contig_match
    elif len(nodes_in_path) == 1:
        gene_contig_match.nodes_list = nodes_in_path
        gene_contig_match.start_in_first_node = contig_start
    else:
        prev_seq = str(records_dict[nodes_in_path[0]].seq)
        end = len(prev_seq)
        if contig_start <= end:  # in case that the match starts in the first node
            start_in_first_node = contig_start
            nodes_for_genes = [nodes_in_path[0]]
        else:
            nodes_for_genes = []
            start_in_first_node = None
        for node in nodes_in_path[1:]:
            cur_seq = str(records_dict[node].seq)
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


def get_genes_to_contigs_with_nodes_list(genes_path, contigs_path, temp_files_path, assembly_graph, records_dict,
                                         paths_path, n_minimap_threads, pident_filtering_th):
    genes_to_contigs_path = f'{temp_files_path}/genes_to_contigs.paf'

    genes_to_contigs_path = sau.map_genes_to_contexts(genes_path, contigs_path, genes_to_contigs_path,
                                                      nthreads=n_minimap_threads)
    genes_to_contigs = sau.read_and_filter_minimap_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                           pident_filtering_th,
                                                           verbose=True)
    contig_nodes_dict = pu.get_contig_nodes_dict(assembly_graph.nodes, paths_path, keep_contigs_with_gaps=False)
    matches_with_nodes_list_and_start_location = []
    for gene_contig_match in genes_to_contigs:  # TODO turn this to an iterator when I finish debugging
        updated_gene_contig_match = add_nodes_list_and_start_location_to_gene_contig_match(records_dict,
                                                                                           contig_nodes_dict.get(
                                                                                               gene_contig_match.contig,
                                                                                               None), gene_contig_match)
        if updated_gene_contig_match.nodes_list and updated_gene_contig_match.start_in_first_node:
            matches_with_nodes_list_and_start_location.append(updated_gene_contig_match)
    return matches_with_nodes_list_and_start_location


def find_genes_in_graph(assembly_graph_path, contigs_path, gene_pident_filtering_th, genes_path, n_minimap_threads,
                        paths_path, temp_folder):
    assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    records_dict = pu.get_records_dict_from_assembly_graph(assembly_graph_path)
    log.info(f'{dt.datetime.now()} mapping genes to contigs')
    genes_to_contigs = get_genes_to_contigs_with_nodes_list(genes_path, contigs_path, temp_folder, assembly_graph,
                                                            records_dict, paths_path, n_minimap_threads,
                                                            gene_pident_filtering_th)
    return assembly_graph, genes_to_contigs, records_dict
