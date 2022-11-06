import os
from Bio import SeqIO
import logging

log = logging.getLogger(__name__)


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


def get_records_dict_from_assembly_graph(assembly_graph_path):
    with open(assembly_graph_path) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    records_dict = {get_short_node_name(record.name): record for record in records}
    return records_dict


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


def is_similar_to_representatives(representatives, record, iou_th):
    for rep in representatives:
        if interval_iou(record.start, record.end, rep.start, rep.end) > iou_th:
            return True
    return False
