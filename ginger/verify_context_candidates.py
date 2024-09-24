import logging
from collections import defaultdict
import datetime as dt

from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu

from typing import Dict, List
log = logging.getLogger(__name__)

# TODO go over all functions and see if they are used and or should be re written
# TODO proper logging
IOU_TH = 0.5


def extract_start_and_end(in_match: mc.PathRefGenomeMatch, out_match: mc.PathRefGenomeMatch) -> tuple:
    if in_match.strand == '+':
        start = in_match.ref_genome_end
        end = out_match.ref_genome_start
    else:
        start = out_match.ref_genome_end
        end = in_match.ref_genome_start
    return start, end


def get_in_out_match(i, o, gene_length, minimal_gap_ratio, maximal_gap_ratio):
    for field in ['gene', 'reference_genome', 'strand']:
        if getattr(i, field) != getattr(o, field):
            return None
    start, end = extract_start_and_end(i, o)
    start_end_diff = end - start
    gap_ratio = start_end_diff / gene_length
    score = (i.score * i.path_length + o.score * o.path_length) / (i.path_length + o.path_length)
    if minimal_gap_ratio < gap_ratio < maximal_gap_ratio:
        return mc.InOutPathsMatch(i, o, start, end, gap_ratio, score, gene_length)
    return None


def read_and_filter_path_matches_per_gene(match_object_constructor: callable, alignment_path, pident_filtering_th, ref_species_dict):
    parsed_as_iterator = sau.read_and_filter_minimap_matches(match_object_constructor, alignment_path,
                                                             pident_filtering_th, ref_species_dict)
    if parsed_as_iterator is None:
        return []
    genes_to_matches = defaultdict(list)
    for match in parsed_as_iterator:
        # print(match)
        genes_to_matches[(match.gene, match.ref_genome)].append(match)
    log.info(
        f"found {sum(len(v) for v in genes_to_matches.values())} matches for {len(genes_to_matches)} gene and ref genomes pairs")
    return genes_to_matches


def get_all_in_out_matches(in_paths_by_gene_and_ref_genome, out_paths_by_gene_and_ref_genome, genes_lengths, minimal_gap_ratio,
                           maximal_gap_ratio, iou_th=IOU_TH) -> Dict[tuple, list]:
    matches_per_gene_and_ref_genome = dict()
    for gene_ref_genome in in_paths_by_gene_and_ref_genome:
        matches_for_gene_ref_genome_pair = []
        in_paths = in_paths_by_gene_and_ref_genome.get(gene_ref_genome, [])
        out_paths = out_paths_by_gene_and_ref_genome.get(gene_ref_genome, [])
        for i in in_paths:
            for o in out_paths:
                in_out_match = get_in_out_match(i, o, genes_lengths[i.gene], minimal_gap_ratio,maximal_gap_ratio)
                if in_out_match is not None:
                    matches_for_gene_ref_genome_pair.append(in_out_match)
        if matches_for_gene_ref_genome_pair:
            matches_per_gene_and_ref_genome[gene_ref_genome]= keep_best_matches(matches_for_gene_ref_genome_pair, iou_th=iou_th)
    log.info(
        f"found {sum((len(m) for m in matches_per_gene_and_ref_genome.values()))} matches for {len(matches_per_gene_and_ref_genome)} gene-reference-genome pairs")
    return matches_per_gene_and_ref_genome

def keep_best_matches(matches:List, sorting_func=lambda x: (-x.score, -(x.end - x.start), x.start),
                      iou_th=IOU_TH):  # -> Dict[Tuple[str,str]:List[mc.InOutPathsMatch]]
    # TODO double check that I can get more than one match
    sorted_matches = sorted(matches, key=sorting_func)
    representative_matches = []
    for match in sorted_matches:
        if not pu.is_similar_to_representatives(representative_matches, match, iou_th):
            representative_matches.append(match)
    return representative_matches


def get_ref_genome_species_dict_from_metadata_path(metadata_path):
    ref_genome_species_dict = {}
    with open(metadata_path, 'r') as f:
        f.readline()
        for line in f:
            splt = line[:-1].split('\t')
            if len(splt) ==4:
                genome, _, _, species = splt
            else:
                print(line)
            ref_genome_species_dict[genome] = species
    return ref_genome_species_dict

@pu.step_timing
def process_in_and_out_paths_to_results(in_path_mapping_to_ref_genomes, out_path_mapping_to_ref_genomes, genes_lengths,
                                        paths_pident_filtering_th, minimal_gap_ratio,
                                        maximal_gap_ratio, metadata_path):
    # TODO get rid of pandas here (the tables have millions of entries and can potentially grow bigger)
    log.info(f'{dt.datetime.now()} parsing the mapping of in and out paths')
    ref_species_dict = get_ref_genome_species_dict_from_metadata_path(metadata_path)
    parsed_in_path_to_ref_genomes_by_gene_and_ref_genome = read_and_filter_path_matches_per_gene(mc.PathRefGenomeMatch,
                                                                                   in_path_mapping_to_ref_genomes,
                                                                                   paths_pident_filtering_th,
                                                                                   ref_species_dict)
    parsed_out_path_to_ref_genomes_by_gene_and_ref_genome = read_and_filter_path_matches_per_gene(mc.PathRefGenomeMatch,
                                                                                    out_path_mapping_to_ref_genomes,
                                                                                    paths_pident_filtering_th,
                                                                                    ref_species_dict)
    if len(parsed_in_path_to_ref_genomes_by_gene_and_ref_genome) == 0 or len(parsed_out_path_to_ref_genomes_by_gene_and_ref_genome) == 0:
        log.info(
            f'GInGeR found {len(parsed_in_path_to_ref_genomes_by_gene_and_ref_genome)=} matches for incoming paths and {len(parsed_out_path_to_ref_genomes_by_gene_and_ref_genome)=} matches for outgoing paths. No results will be produced')
        return []
    log.info(f'{dt.datetime.now()} generating in-out matches')
    matches_per_gene_and_ref_genome = get_all_in_out_matches(parsed_in_path_to_ref_genomes_by_gene_and_ref_genome,
                                              parsed_out_path_to_ref_genomes_by_gene_and_ref_genome,
                                              genes_lengths,
                                              minimal_gap_ratio, maximal_gap_ratio)
    return matches_per_gene_and_ref_genome
