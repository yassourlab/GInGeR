import logging
from collections import defaultdict
import datetime as dt

from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu

log = logging.getLogger(__name__)


# TODO go over all functions and see if they are used and or should be re written
# TODO proper logging

def extract_start_and_end(i_start, i_end, o_start, o_end, strand):
    if strand == '+':
        start = i_end
        end = o_start
    else:
        start = o_end
        end = i_start
    return start, end


def get_in_out_match(i, o, gene_length, minimal_gap_ratio, maximal_gap_ratio):
    for field in ['gene', 'nodes_list', 'bug', 'strand']:
        if getattr(i, field) != getattr(o, field):
            return None
    start, end = extract_start_and_end(i.bug_start, i.bug_end, o.bug_start, o.bug_end, i.strand)
    start_end_diff = end - start
    gap_ratio = start_end_diff / gene_length
    match_score = (i.match_score * i.path_length + o.match_score * o.path_length) / (i.path_length + o.path_length)
    if minimal_gap_ratio < gap_ratio < maximal_gap_ratio:
        return mc.InOutPathsMatch(i, o, start, end, gap_ratio, match_score, gene_length)
    return None


def read_and_filter_path_matches_per_gene(match_object_constructor: callable, alignment_path, pident_filtering_th):
    parsed_as_iterator = sau.read_and_filter_minimap_matches(match_object_constructor, alignment_path,
                                                             pident_filtering_th)
    if parsed_as_iterator is None:
        return []
    genes_to_matches = defaultdict(list)
    for match in parsed_as_iterator:
        genes_to_matches[(match.gene, match.bug)].append(match)
    log.info(
        f"found {sum(len(v) for v in genes_to_matches.values())} matches for {len(genes_to_matches)} gene-bug pairs")
    return genes_to_matches


def get_all_in_out_matches(in_paths_by_gene_and_bug, out_paths_by_gene_and_bug, genes_lengths, minimal_gap_ratio,
                           maximal_gap_ratio) -> defaultdict:
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
    log.info(
        f"found {sum((len(m) for m in matches_per_gene_and_bug.values()))} matches for {len(matches_per_gene_and_bug)} gene-bug pairs")
    return matches_per_gene_and_bug


def keep_best_matches_per_gene_bug_pair(matches_per_gene) -> defaultdict:
    best_matches_per_gene = defaultdict()
    for gene, matches in matches_per_gene.items():
        best_matches_per_gene[gene] = keep_best_matches(matches, iou_th=0.5)
    return best_matches_per_gene


def keep_best_matches(matches, iou_th: float = 0):  # -> Dict[Tuple[str,str]:List[mc.InOutPathsMatch]]
    # TODO double check that I can get more than one match
    sorted_matches = sorted(matches, key=lambda x: x.match_score, reverse=True)
    representative_matches = []
    for match in sorted_matches:
        if not pu.is_similar_to_representatives(representative_matches, match, iou_th):
            representative_matches.append(match)
    return representative_matches


def process_in_and_out_paths_to_results(in_path_mapping_to_bugs, out_path_mapping_to_bugs, genes_lengths,
                                        paths_pident_filtering_th, minimal_gap_ratio,
                                        maximal_gap_ratio):
    # TODO get rid of pandas here (the tables have millions of entries and can potentially grow bigger)
    log.info(f'{dt.datetime.now()} parsing the mapping of in and out paths')
    parsed_in_path_to_bugs_by_gene_and_bug = read_and_filter_path_matches_per_gene(mc.PathBugMatch,
                                                                                   in_path_mapping_to_bugs,
                                                                                   paths_pident_filtering_th)
    parsed_out_path_to_bugs_by_gene_and_bug = read_and_filter_path_matches_per_gene(mc.PathBugMatch,
                                                                                    out_path_mapping_to_bugs,
                                                                                    paths_pident_filtering_th)
    if len(parsed_in_path_to_bugs_by_gene_and_bug)==0 or len(parsed_out_path_to_bugs_by_gene_and_bug)==0:
        log.info(f'GInGeR found {len(parsed_in_path_to_bugs_by_gene_and_bug)=} matches for incoming paths and {len(parsed_out_path_to_bugs_by_gene_and_bug)=} matches for outgoing paths. No results will be produced')
        return []
    log.info(f'{dt.datetime.now()} generating in-out matches')
    matches_per_gene = get_all_in_out_matches(parsed_in_path_to_bugs_by_gene_and_bug,
                                              parsed_out_path_to_bugs_by_gene_and_bug,
                                              genes_lengths,
                                              minimal_gap_ratio, maximal_gap_ratio)
    log.info(f'{dt.datetime.now()} filtering in-out matches')
    matches_per_gene_no_overlaps = keep_best_matches_per_gene_bug_pair(matches_per_gene)
    return matches_per_gene_no_overlaps
