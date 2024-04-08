import os
import logging
import pandas as pd
import timeit
from collections import defaultdict

RUNTIME_PRINTS_PATTERN = '$$$$$$$$$$'
log = logging.getLogger(__name__)

# TODO - write a function here that runs an external tool and present the output using tqdm (I coppied and pasted it multiple times already)
def step_timing(func):
    def wrapper_lot_and_time(*args):
        start = timeit.default_timer()
        func_return_vals = func(*args)
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


def parse_paths_file(paths_file, assembly_graph_nodes, path_is_contig_name_func=is_contig_name_func,
                     keep_contigs_with_gaps=True):
    out = []
    contig_name = None
    nodes_details_list = []
    had_gaps = False
    for line in paths_file.readlines():
        line_no_newline = line[:-len('\n')]
        if path_is_contig_name_func(line_no_newline):
            if contig_name is not None and (keep_contigs_with_gaps or not had_gaps):
                out.append((contig_name, nodes_details_list))
                # yield contig_name, nodes_details_list
            contig_name = line_no_newline
            nodes_details_list = []
            had_gaps = False
        else:
            had_gaps = had_gaps or (';' in line)
            nodes_details_list += parse_list_of_nodes(line_no_newline, assembly_graph_nodes)
    out.append((contig_name, nodes_details_list))
    return out


def get_contig_nodes_dict(assembly_graph_nodes, paths_path, keep_contigs_with_gaps=True):
    with open(paths_path) as f:
        # TODO - I just ignore contigs with gaps and I can technically miss matches this way
        contig_nodes_iterator = parse_paths_file(f, assembly_graph_nodes, keep_contigs_with_gaps=keep_contigs_with_gaps)
    contig_nodes_dict = {contig_name: nodes_list for contig_name, nodes_list in contig_nodes_iterator}
    return contig_nodes_dict


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


def interval_iou(start1, end1, start2, end2):
    try:
        # TODO can probably calculate intersection in a shorter way (see tal's snippet at WW) +  move it elsewhere
        if start1 <= start2 <= end1 <= end2:

            intersection = end1 - start2
            union = end2 - start1
        elif start2 <= start1 <= end2 <= end1:

            intersection = end2 - start1
            union = end1 - start1

        elif start1 <= start2 <= end2 <= end1:

            intersection = end2 - start2
            union = end1 - start1
        elif start2 <= start1 <= end1 <= end2:

            intersection = end1 - start1
            union = end2 - start2
        else:
            return 0
        output = intersection / union
    except Exception as e:
        log.error(f'{e} - start1: {start1}, end1: {end1}, start2: {start2}, end2: {end2}')
        raise e
    return output


def is_similar_to_representatives(representatives, gene_paths_to_bug_match, iou_th):
    for rep in representatives:
        if interval_iou(gene_paths_to_bug_match.start, gene_paths_to_bug_match.end, rep.start, rep.end) > iou_th:
            return True
    return False

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
                results_dict['in_context_start'].append(match.in_path.bug_start)
                results_dict['out_context_end'].append(match.out_path.bug_end)
            else:
                results_dict['in_context_start'].append(match.out_path.bug_start)
                results_dict['out_context_end'].append(match.in_path.bug_end)
            results_dict['gene_start'].append(match.start)
            results_dict['gene_end'].append(match.end)

            results_dict['score'].append(match.score)

    metadata_df = pd.read_csv(metadata_path, sep='\t')
    results_df = pd.DataFrame(results_dict)
    results_df['Genome'] = results_df['reference_contig'].apply(lambda x: x.split('_')[0].split('.')[0])
    results_df.merge(metadata_df[['Genome', 'species']], on='Genome', how='left').to_csv(csv_path, index=False)

@step_timing
def aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                         species_level_output_path,
                                                                         max_species_representatives: int):
    context_level_df = pd.read_csv(context_level_output_path)
    metadata_df = pd.read_csv(metadata_path, sep='\t')
    context_level_df_with_metadata = pd.merge(context_level_df, metadata_df[['Genome', 'species']],
                                              on=['Genome', 'species'],
                                              how='left')

    genomes_per_species = metadata_df.species.value_counts().to_frame()
    genomes_per_species['species_instances'] = genomes_per_species['count'].apply(
        lambda x: min(x, max_species_representatives))
    agg_output = context_level_df_with_metadata.groupby(['gene', 'species']).aggregate(
        {'Genome': ['nunique'], 'score': ['max']})
    agg_output.columns = ['_'.join(col) for col in agg_output.columns.values]
    agg_output = agg_output.merge(genomes_per_species, left_on='species', right_index=True, how='left')
    agg_output['references_ratio'] = agg_output['Genome_nunique'] / agg_output['species_instances']
    agg_output[['references_ratio', 'score_max', 'species_instances']].to_csv(species_level_output_path)
