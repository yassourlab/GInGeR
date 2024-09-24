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


def is_similar_to_representatives(representatives, gene_paths_to_ref_genome_match, iou_th):
    for rep in representatives:
        if interval_iou(gene_paths_to_ref_genome_match.start, gene_paths_to_ref_genome_match.end, rep.start, rep.end) > iou_th:
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
                results_dict['in_context_start'].append(match.in_path.ref_genome_start)
                results_dict['out_context_end'].append(match.out_path.ref_genome_end)
            else:
                results_dict['in_context_start'].append(match.out_path.ref_genome_start)
                results_dict['out_context_end'].append(match.in_path.ref_genome_end)
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

    context_level_df.columns = [x.lower() for x in context_level_df.columns]
    metadata_df.columns = [x.lower() for x in metadata_df.columns]

    context_level_df_with_metadata = pd.merge(context_level_df, metadata_df[['genome', 'species']],
                                              on=['genome', 'species'],
                                              how='left')

    genomes_per_species = metadata_df.species.value_counts().to_frame()
    genomes_per_species['species_instances'] = genomes_per_species['count'].apply(
        lambda x: min(x, max_species_representatives))
    agg_output = context_level_df_with_metadata.groupby(['gene', 'species']).aggregate(
        {'genome': ['nunique'], 'score': ['max']})
    agg_output.columns = ['_'.join(col) for col in agg_output.columns.values]
    agg_output = agg_output.merge(genomes_per_species, left_on='species', right_index=True, how='left')
    agg_output['references_ratio'] = agg_output['genome_nunique'] / agg_output['species_instances']
    species_level_output = agg_output[['references_ratio', 'score_max', 'species_instances']]
    if species_level_output_path is not None:
        species_level_output.to_csv(species_level_output_path)
    return species_level_output
