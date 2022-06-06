import math

# import pandas as pd
import pyfastg
import datetime as dt
import networkx as nx

import context_graph_approach as cga
import sequence_alignment_utils as sau
import contigs_alignment_approach as caa
from collections import defaultdict

# TODO go over all functions and see if they are used and or should be re written
# TODO proper logging
# # using pandas - keeping this for generating UHGG reference
# def extract_start_and_end_series(i_start, i_end, o_start, o_end, strand):
#     if strand == '+':
#         return pd.Series({'start': i_end, 'end': o_start})
#     else:
#         return pd.Series({'start': o_end, 'end': i_start})


# def extract_results_from_in_and_out_mappings(all_in_out_paths, minimal_gap_ratio, maximal_gap_ratio,
#                                              temp_file_path=None, drop_dups=True):
#     if len(all_in_out_paths) == 0:
#         return pd.DataFrame()
#
#     start_end_df = all_in_out_paths.apply(
#         lambda x: extract_start_and_end_series(x.bug_start_i, x.bug_end_i, x.bug_start_o, x.bug_end_o, x.strand),
#         axis=1)
#     all_in_out_paths = all_in_out_paths.merge(start_end_df, left_index=True, right_index=True)
#     all_in_out_paths['start_end_diff'] = all_in_out_paths.end - all_in_out_paths.start
#
#     all_in_out_paths['gap_ratio'] = (all_in_out_paths.start_end_diff / all_in_out_paths.gene_length)
#     # all_in_out_paths['gap_ratio'] = (all_in_out_paths.start_end_diff / all_in_out_paths.gene_nodes_length)
#     path_pident_numerator = all_in_out_paths.matches_ratio_i * all_in_out_paths.contig_length_i + all_in_out_paths.matches_ratio_o * all_in_out_paths.contig_length_o
#     path_pident_denomenator = all_in_out_paths.contig_length_i + all_in_out_paths.contig_length_o
#     all_in_out_paths['match_score'] = path_pident_numerator / path_pident_denomenator
#     if temp_file_path is not None:
#         all_in_out_paths.to_csv(f'{temp_file_path}/all_in_out_paths.csv',
#                                 index=False)
#     passed_filtering = all_in_out_paths[
#         (all_in_out_paths.gap_ratio > minimal_gap_ratio) & (
#                 all_in_out_paths.gap_ratio < maximal_gap_ratio)].sort_values(
#         'match_score')
#     if drop_dups:  # TODO I'm not sure why I did this dup removal in the first place...
#         passed_filtering = passed_filtering.drop_duplicates(
#             ['bug', 'gene', 'bug_start_i'], keep='last').drop_duplicates(
#             ['bug', 'gene', 'bug_start_o'], keep='last').drop_duplicates(
#             ['bug', 'gene', 'bug_end_i'], keep='last').drop_duplicates(
#             ['bug', 'gene', 'bug_end_o'], keep='last')
#
#     return passed_filtering


# The up to date code for this approach
def extract_start_and_end(i_start, i_end, o_start, o_end, strand):
    if strand == '+':
        start = i_end
        end = o_start
    else:
        start = o_end
        end = i_start
    return start, end


def save_paths_to_fasta_io_paths_approach(paths, paths_fasta_name, records_dict, node_to_find=None,
                                          context_len=math.inf, in_or_out=None, covered_by_gene=0, gene_and_node=''):
    # TODO I think I don't really need the outputs of this function anymore
    # print(f'{dt.datetime.now()} saving paths to {paths_fasta_name}')
    in_paths_lengths = {}
    node_locations = {}
    with open(paths_fasta_name, 'a') as f:
        for path in paths:
            seq, node_start = cga.generate_str_from_list_of_nodes(records_dict, path, node_to_find)
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
                k = cga.get_sequence_overlap(prev_seq, cur_seq)
            except Exception as e:
                print(f'error! {str(e)} {node} ')
                raise Exception
            start = end - k
            end = start + len(cur_seq)
            if cga.intervals_overlap(start, end, contig_start, contig_end):
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
    # processed_genes_to_contigs_path = f'{temp_files_path}/processed_genes_to_contigs.csv'

    genes_to_contigs_path = sau.map_genes_to_contigs(genes_path, contigs_path, genes_to_contigs_path,
                                                     nthreads=n_minimap_threads)
    genes_to_contigs = sau.read_and_filter_minimap_matches(sau.GeneContigMatch, genes_to_contigs_path,
                                                           pident_filtering_th,
                                                           verbose=True)
    contig_nodes_dict = cga.get_contig_nodes_dict(assembly_graph.nodes, paths_path, keep_contigs_with_gaps=False)
    matches_with_nodes_list_and_start_location = []
    for gene_contig_match in genes_to_contigs:  # TODO turn this to an iterator when I finish debugging
        updated_gene_contig_match = add_nodes_list_and_start_location_to_gene_contig_match(records_dict,
                                                                                           contig_nodes_dict.get(
                                                                                               gene_contig_match.contig,
                                                                                               None), gene_contig_match)
        if updated_gene_contig_match.nodes_list and updated_gene_contig_match.start_in_first_node:
            matches_with_nodes_list_and_start_location.append(updated_gene_contig_match)
    print(len(matches_with_nodes_list_and_start_location))
    return matches_with_nodes_list_and_start_location


def extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, records_dict, genes_to_contigs, depth_limit,
                                                      context_len,
                                                      in_paths_fasta, out_paths_fasta):
    gene_and_nodes_path_set = set()
    gene_lengths = {}
    for gene_contigs_match in genes_to_contigs:
        gene_and_nodes_path_str = f"{gene_contigs_match.gene}_nodes_{'_'.join(gene_contigs_match.nodes_list)}"
        if gene_and_nodes_path_str in gene_and_nodes_path_set or gene_contigs_match.start_in_first_node is None:  # I already ran the pipeline for this and there is no need to do it again
            print(
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
            gene_nodes_length = len(cga.generate_str_from_list_of_nodes(records_dict, nodes_list_for_gene, None)[0])

            # in paths
            in_covered_by_gene = assembly_graph.nodes[first_node]['length'] - start_in_first_node
            in_paths_initial_stack = [(first_node, [(first_node, assembly_graph.nodes[first_node]['length'])])]
            in_paths = paths_enumerator(assembly_graph, in_paths_initial_stack, depth_limit, context_len,
                                        nx.DiGraph.predecessors, reverse=True, covered_by_gene=in_covered_by_gene)
            in_paths_lengths, _ = save_paths_to_fasta_io_paths_approach(in_paths, in_paths_fasta, records_dict,
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
            out_paths_lengths, _ = save_paths_to_fasta_io_paths_approach(out_paths, out_paths_fasta, records_dict,
                                                                         context_len=context_len, in_or_out='out',
                                                                         covered_by_gene=out_covered_by_gene,
                                                                         gene_and_node=gene_and_nodes_path_str)
            gene_lengths[gene_name] = gene_length
            if len(in_paths) == 0 or len(out_paths) == 0:
                print(f'{gene_contigs_match} in {len(in_paths)} out {len(out_paths)}')

    return gene_lengths


def run_in_out_paths_approach(assembly_dir, genes_path, reference_path,
                              n_minimap_threads, depth_limit, maximal_gap_ratio, context_len,
                              gene_pident_filtering_th,
                              paths_pident_filtering_th, temp_folder):
    contigs_path = f"{assembly_dir}/{'contigs.fasta'}"
    paths_path = f"{assembly_dir}/{'contigs.paths'}"
    assembly_graph_path = f"{assembly_dir}/{'assembly_graph.fastg'}"
    in_paths_fasta = f"{temp_folder}/all_in_paths.fasta"
    out_paths_fasta = f"{temp_folder}/all_out_paths.fasta"
    in_mapping_to_bugs_path = f'{temp_folder}/in_paths_to_reference.paf'
    out_mapping_to_bugs_path = f'{temp_folder}/out_paths_to_reference.paf'

    assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    records_dict = cga.get_records_dict_from_assembly_graph(assembly_graph_path)
    print(f'{dt.datetime.now()} mapping genes to contigs')
    genes_to_contigs = get_genes_to_contigs_with_nodes_list(genes_path, contigs_path, temp_folder, assembly_graph,
                                                            records_dict, paths_path, n_minimap_threads,
                                                            gene_pident_filtering_th)
    # TODO if there are no in paths, there is no need to calc out paths. but maybe we should always have in paths. check this
    # get in and out paths
    print(f'{dt.datetime.now()} extracting in and out paths')
    genes_lengths = extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, records_dict,
                                                                      genes_to_contigs, depth_limit,
                                                                      context_len,
                                                                      in_paths_fasta, out_paths_fasta)
    # map them to the reference
    print(f'{dt.datetime.now()} mapping in and out paths')
    sau.map_contigs_to_bugs(in_paths_fasta, reference_path, in_mapping_to_bugs_path, nthreads=n_minimap_threads)
    sau.map_contigs_to_bugs(out_paths_fasta, reference_path, out_mapping_to_bugs_path, nthreads=n_minimap_threads)

    # merge and get results
    return process_in_and_out_paths_to_results_no_pd(in_mapping_to_bugs_path, out_mapping_to_bugs_path, genes_lengths,
                                                     paths_pident_filtering_th, 0, maximal_gap_ratio)


class InOutPathsMatch:
    def __init__(self, in_path, out_path, start, end, gap_ratio, match_score, gene_length, gene=None, bug=None):
        self.in_path = in_path
        self.out_path = out_path
        if gene is None:
            self.gene = self.in_path.gene
        else:
            self.gene = gene

        if bug is None:
            self.bug = self.in_path.bug
        else:
            self.bug = bug

        self.start = start
        self.end = end
        self.gap_ratio = gap_ratio
        self.match_score = match_score
        self.gene_length = gene_length

    def to_dict(self):
        out_dict = dict(gene=self.gene, bug=self.bug, start=self.start, end=self.end, match_score=self.match_score)
        if self.in_path:
            out_dict.update(dict(in_path_name=self.in_path.query_name, out_path_name=self.out_path.query_name,
                                 in_path_score=self.in_path.match_score, out_path_score=self.out_path.match_score))

        return out_dict


def get_in_out_match(i, o, gene_length, minimal_gap_ratio, maximal_gap_ratio):
    for field in ['gene', 'nodes_list', 'bug', 'strand']:
        if getattr(i, field) != getattr(o, field):
            return None
    start, end = extract_start_and_end(i.bug_start, i.bug_end, o.bug_start, o.bug_end, i.strand)
    start_end_diff = end - start
    gap_ratio = start_end_diff / gene_length
    match_score = (i.match_score * i.path_length + o.match_score * o.path_length) / (i.path_length + o.path_length)
    if minimal_gap_ratio < gap_ratio < maximal_gap_ratio:
        return InOutPathsMatch(i, o, start, end, gap_ratio, match_score, gene_length)
    return None


def read_and_filter_path_matches_per_gene(match_object_constructor: callable, alignment_path, pident_filtering_th,
                                          verbose=True):
    parsed_as_iterator = sau.read_and_filter_minimap_matches(match_object_constructor, alignment_path,
                                                             pident_filtering_th,
                                                             False)
    genes_to_matches = defaultdict(list)
    for match in parsed_as_iterator:
        genes_to_matches[(match.gene, match.bug)].append(match)
    if verbose:
        print(
            f"found {sum(len(v) for v in genes_to_matches.values())} matches for {len(genes_to_matches)} gene-bug pairs")
    return genes_to_matches


def get_all_in_out_matches(in_paths_by_gene_and_bug, out_paths_by_gene_and_bug, genes_lengths, minimal_gap_ratio,
                           maximal_gap_ratio):
    matches_per_gene_and_bug = defaultdict(list)
    for gene_bug in in_paths_by_gene_and_bug:
        if gene_bug in out_paths_by_gene_and_bug:
            in_paths = in_paths_by_gene_and_bug[gene_bug]
            out_paths = out_paths_by_gene_and_bug[gene_bug]
            for i in in_paths:
                for o in out_paths:
                    in_out_match = get_in_out_match(i, o, genes_lengths[i.gene], minimal_gap_ratio, maximal_gap_ratio)
                    if in_out_match is not None:
                        matches_per_gene_and_bug[gene_bug].append(in_out_match)
    print(
        f"found {sum((len(m) for m in matches_per_gene_and_bug.values()))} matches for {len(matches_per_gene_and_bug)} gene-bug pairs")
    return matches_per_gene_and_bug


def keep_best_matches_per_gene_bug_pair(matches_per_gene):
    best_matches_per_gene = defaultdict()
    for gene, matches in matches_per_gene.items():
        best_matches_per_gene[gene] = keep_best_matches(matches, iou_th=0.5)
    return best_matches_per_gene


def keep_best_matches(matches, iou_th: float = 0):
    sorted_matches = sorted(matches, key=lambda x: x.match_score, reverse=True)
    representative_matches = []
    for match in sorted_matches:
        if not caa.is_similar_to_representatives(representative_matches, match, iou_th):
            representative_matches.append(match)
    return representative_matches


def process_in_and_out_paths_to_results_no_pd(in_path_mapping_to_bugs, out_path_mapping_to_bugs, genes_lengths,
                                              paths_pident_filtering_th, minimal_gap_ratio,
                                              maximal_gap_ratio):
    # TODO get rid of pandas here (the tables have millions of entries and can potentially grow bigger)
    print(f'{dt.datetime.now()} parsing the mapping of in and out paths')
    parsed_in_path_to_bugs_by_gene_and_bug = read_and_filter_path_matches_per_gene(sau.PathBugMatch,
                                                                                   in_path_mapping_to_bugs,
                                                                                   paths_pident_filtering_th,
                                                                                   verbose=True)
    parsed_out_path_to_bugs_by_gene_and_bug = read_and_filter_path_matches_per_gene(sau.PathBugMatch,
                                                                                    out_path_mapping_to_bugs,
                                                                                    paths_pident_filtering_th,
                                                                                    verbose=True)
    print(f'{dt.datetime.now()} generating in-out matches')
    matches_per_gene = get_all_in_out_matches(parsed_in_path_to_bugs_by_gene_and_bug,
                                              parsed_out_path_to_bugs_by_gene_and_bug,
                                              genes_lengths,
                                              minimal_gap_ratio, maximal_gap_ratio)
    print(f'{dt.datetime.now()} filtering in-out matches')
    matches_per_gene_no_overlaps = keep_best_matches_per_gene_bug_pair(matches_per_gene)
    return matches_per_gene_no_overlaps

