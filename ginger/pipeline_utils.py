import os
import logging
import pandas as pd
import timeit
from collections import defaultdict
from pafpy import PafRecord, PafFile

RUNTIME_PRINTS_PATTERN = '$$$$$$$$$$'
log = logging.getLogger(__name__)


# TODO - write a function here that runs an external tool and present the output using tqdm (I coppied and pasted it multiple times already)
def step_timing(func):
    def wrapper_lot_and_time(*args, **kwargs):
        start = timeit.default_timer()
        func_return_vals = func(*args, **kwargs)
        stop = timeit.default_timer()
        log.info(f'{RUNTIME_PRINTS_PATTERN} {func} took {(stop - start) / 60} minutes {RUNTIME_PRINTS_PATTERN}')
        return func_return_vals

    return wrapper_lot_and_time


def check_and_makedir(path_with_file):
    path = '/'.join(path_with_file.split('/')[:-1])
    if not os.path.exists(path):
        os.makedirs(path)


def check_and_make_dir_no_file_name(path):
    if not os.path.exists(path):
        os.makedirs(path)


def generate_str_from_list_of_nodes(records_dict, nodes_in_path, node_to_find=None):
    node_start, node_end = None, None
    prev_seq = str(records_dict[nodes_in_path[0]].seq)
    prev_node = nodes_in_path[0]
    seq = prev_seq[:]
    for node in nodes_in_path[1:]:
        cur_seq = str(records_dict[node].seq)
        try:
            k = get_sequence_overlap(prev_seq, cur_seq)
        except Exception as e:
            logging.error(f'{str(e)} {prev_node} {node} ')
            raise Exception
        if node == node_to_find:  # assuming that the node appears only once. will take the last location of the node
            node_start = len(seq) - k

        seq += cur_seq[k:]
        prev_node = node
        prev_seq = cur_seq
    return seq, node_start


def parse_list_of_nodes(as_str, original_graph_nodes):
    splt = as_str.split(',')
    return [node.replace(';', '') for node in splt]


def is_contig_name_func(line):
    return line.startswith('NODE')


def paf_record_to_dict(paf_record):
    return dict(qname=paf_record.qname, qlen=paf_record.qlen, strand=str(paf_record.strand),
                tname=paf_record.tname, tlen=paf_record.tlen, tstart=paf_record.tstart, tend=paf_record.tend,
                mlen=paf_record.mlen)  # blen=paf_record.blen, qstart=paf_record.qstart, qend=paf_record.qend,


def minimap_results_from_path(path, head_size=None):
    with open(path) as f:
        paf_file = PafFile(f)
        if head_size:
            head = [next(paf_file) for _ in range(head_size)]  # paf_file
        else:
            head = paf_file
        minimap_results = pd.DataFrame([paf_record_to_dict(paf_record) for paf_record in head])
    return minimap_results


def parse_paths_file(paths_path, assembly_graph_nodes, path_is_contig_name_func=is_contig_name_func):
    parsed_paths_without_gaps = {}
    contigs_with_gaps = set()

    contig_name = None
    path_in_graph = []
    had_gaps = False
    with open(paths_path) as paths_file:
        for line in paths_file.readlines():
            line_no_newline = line[:-len('\n')]
            if path_is_contig_name_func(line_no_newline):
                if contig_name:  # add contig name to the dict of parsed paths or to the list of paths with gaps
                    if had_gaps:
                        contigs_with_gaps.add(contig_name)
                    else:
                        parsed_paths_without_gaps[contig_name] = path_in_graph

                contig_name = line_no_newline
                path_in_graph = []
                had_gaps = False
            else:
                had_gaps = had_gaps or (';' in line)
                if not had_gaps:
                    path_in_graph += parse_list_of_nodes(line_no_newline, assembly_graph_nodes)
        # add the last contig
        if had_gaps:
            contigs_with_gaps.add(contig_name)
        else:
            parsed_paths_without_gaps[contig_name] = path_in_graph
        return parsed_paths_without_gaps, contigs_with_gaps


def get_sequence_overlap(seq_a, seq_b):
    ks = [55, 33, 21, 43]  # I added 43 because I found it myself. it didn't appear in the documentation
    for k in ks:
        if seq_a[-k:] == seq_b[:k]:
            return k
    possible_k = None
    for k in range(min([len(seq_a), len(seq_b), 200]), 0, -1):
        if seq_a[-k:] == seq_b[:k]:
            possible_k = k
            break
    raise Exception(
        f'No overlap was found with ks {ks}. possible k between 1 to {min([len(seq_a), len(seq_b), 200])} - k={possible_k}. {seq_a}\n {seq_b}')


def intervals_overlap(start_a, end_a, start_b, end_b):
    return start_a <= end_b and start_b <= end_a


def write_genes_detected_in_graph_with_no_species_match(genes_with_location_in_graph, matched_genes, csv_path: str) -> bool:
    """Write a CSV listing detected genes in the assembly graph that lack a species-level match.

    The output has exactly these columns:
    - gene
    - contig
    - gene_match_score

    The file is written only if at least one unmatched gene exists.

    Returns True if the file was written, False otherwise.
    """
    if not genes_with_location_in_graph:
        return False

    matched_genes = set(matched_genes or [])
    rows = []
    for gene_match in genes_with_location_in_graph:
        if gene_match.gene not in matched_genes:
            rows.append(
                {
                    'gene': gene_match.gene,
                    'contig': gene_match.contig,
                    'gene_match_score': gene_match.score,
                }
            )

    if not rows:
        return False

    pd.DataFrame(rows, columns=['gene', 'contig', 'gene_match_score']).to_csv(csv_path, index=False)
    return True





@step_timing
def write_context_level_output_to_csv(output, csv_path: str, metadata_path: str):
    results_dict = defaultdict(list)
    for gene_species_tuple, matches_list in output.items():
        gene, reference = gene_species_tuple
        for match in matches_list:
            results_dict['gene'].append(gene)
            results_dict['reference_contig'].append(reference)
            results_dict['in_context'].append(match.in_path.query_name)
            results_dict['out_context'].append(match.out_path.query_name)
            if match.in_path.strand == '+':
                results_dict['in_context_start'].append(match.in_path.ref_genome_start)
                results_dict['out_context_end'].append(match.out_path.ref_genome_end)
            else:
                results_dict['in_context_start'].append(match.out_path.ref_genome_start)
                results_dict['out_context_end'].append(match.in_path.ref_genome_end)
            results_dict['gene_start'].append(match.start)
            results_dict['gene_end'].append(match.end)

            results_dict['score'].append(match.score)
            results_dict['gene_match_score'].append(match.gene_match_score)
            results_dict['in_context_score'].append(match.in_context_score)
            results_dict['out_context_score'].append(match.out_context_score)

    metadata_df = pd.read_csv(metadata_path, sep='\t')
    results_df = pd.DataFrame(results_dict)
    results_df['Genome'] = results_df['reference_contig'].apply(lambda x: x.split('_')[0].split('.')[0])
    metadata_cols_to_merge = [x for x in ['Genome', 'species','subspecies'] if x in metadata_df.columns]
    results_df.merge(metadata_df[metadata_cols_to_merge], on='Genome', how='left').to_csv(csv_path, index=False)


@step_timing
def aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                         species_level_output_path,
                                                                         max_species_representatives: int,
                                                                         species_col='species'):
    context_level_df = pd.read_csv(context_level_output_path)
    metadata_df = pd.read_csv(metadata_path, sep='\t')

    context_level_df.columns = [x.lower() for x in context_level_df.columns]
    metadata_df.columns = [x.lower() for x in metadata_df.columns]

    genomes_per_species = metadata_df[species_col].value_counts().to_frame()
    genomes_per_species[f'{species_col}_instances'] = genomes_per_species['count'].apply(
        lambda x: min(x, max_species_representatives))
    agg_output = context_level_df.groupby(['gene', species_col]).aggregate(
        {'genome': ['nunique'], 'score': ['max']})
    agg_output.columns = ['_'.join(col) for col in agg_output.columns.values]
    agg_output = agg_output.merge(genomes_per_species, left_on=species_col, right_index=True, how='left')
    agg_output['references_ratio'] = agg_output['genome_nunique'] / agg_output[f'{species_col}_instances']
    species_level_output = agg_output[['references_ratio', 'score_max', f'{species_col}_instances']]

    if 'plasmid_score' in context_level_df.columns:
        group_cols = ['gene', species_col]
        context_cols = ['in_context', 'out_context']

        # average plasmid score across the gene's unique contexts (don't over-weight contexts
        # that were matched to many reference genomes of the same species)
        unique_contexts = context_level_df.drop_duplicates(group_cols + context_cols)
        plasmid_score_mean = unique_contexts.groupby(group_cols)['plasmid_score'].mean().rename('plasmid_score_mean')

        # plasmid score of the context(s) matched to the most reference genomes, averaging ties
        context_counts = context_level_df.groupby(group_cols + context_cols).size().reset_index(name='n_genomes')
        context_counts = context_counts.merge(unique_contexts[group_cols + context_cols + ['plasmid_score']],
                                               on=group_cols + context_cols)
        max_counts = context_counts.groupby(group_cols)['n_genomes'].transform('max')
        plasmid_score_most_common_context = context_counts[context_counts['n_genomes'] == max_counts].groupby(
            group_cols)['plasmid_score'].mean().rename('plasmid_score_most_common_context')

        species_level_output = species_level_output.join(plasmid_score_mean).join(plasmid_score_most_common_context)

    if species_level_output_path is not None:
        species_level_output.to_csv(species_level_output_path)
    return species_level_output
