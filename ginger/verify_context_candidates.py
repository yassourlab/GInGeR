import csv
import logging
import os
from collections import defaultdict
import re

from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu

from typing import Dict, List
log = logging.getLogger(__name__)
# how much two in-out matches on the same reference genome may overlap before the lower scoring one is
# dropped as a redundant description of the same locus. Higher than the threshold gene alignments are
# deduplicated with (sequence_alignment_utils.IOU_TH), because these span a whole context-gene-context
IN_OUT_MATCH_IOU_TH = 0.5


def extract_start_and_end(in_match: mc.PathRefGenomeMatch, out_match: mc.PathRefGenomeMatch) -> tuple:
    if in_match.strand == '+':
        start = in_match.ref_genome_end
        end = out_match.ref_genome_start
    else:
        start = out_match.ref_genome_end
        end = in_match.ref_genome_start
    return start, end


def get_in_out_match(i, o, gene_length, minimal_gap_ratio, maximal_gap_ratio):
    for field in ['gene', 'ref_genome', 'strand', 'locus']:
        if getattr(i, field) != getattr(o, field):
            return None
    start, end = extract_start_and_end(i, o)
    start_end_diff = end - start
    gap_ratio = start_end_diff / gene_length
    score = (i.score * i.path_length + o.score * o.path_length) / (i.path_length + o.path_length)
    if minimal_gap_ratio < gap_ratio < maximal_gap_ratio:
        return mc.InOutPathsMatch(i, o, start, end, gap_ratio, score, gene_length, gene_match_score=i.gene_match_score, in_context_score=i.score, out_context_score=o.score, locus=i.locus)
    return None


def read_and_filter_path_matches_per_gene(match_object_constructor: callable, alignment_path, pident_filtering_th,
                                          ref_species_dict):
    """Groups a context fasta's matches to the reference by the gene copy the context was cut from,
    as well as by the gene and the reference genome."""
    parsed_as_iterator = sau.read_and_filter_minimap_matches(match_object_constructor, alignment_path,
                                                             pident_filtering_th, ref_species_dict)
    if parsed_as_iterator is None:
        return []
    genes_to_matches = defaultdict(list)
    for match in parsed_as_iterator:
        genes_to_matches[(match.gene, match.locus, match.ref_genome)].append(match)
    log.info(
        f"found {sum(len(v) for v in genes_to_matches.values())} matches for {len(genes_to_matches)} gene, gene copy and ref genome triples")
    return genes_to_matches


def get_all_in_out_matches(in_paths_by_gene_locus_and_ref_genome, out_paths_by_gene_locus_and_ref_genome, genes_lengths,
                           minimal_gap_ratio, maximal_gap_ratio, iou_th=IN_OUT_MATCH_IOU_TH) -> Dict[tuple, list]:
    """Pairs every incoming context of a gene copy with every outgoing context of the same copy.

    Pairing across copies would describe a stretch of sequence that is in no contig: the flank of
    one copy, then a gene, then the flank of another. The result is still keyed by (gene, reference
    genome) and deduplicated there, so two copies landing on the same place in the reference are
    still reported once.
    """
    matches_per_gene_and_ref_genome = defaultdict(list)
    for gene, locus, ref_genome in in_paths_by_gene_locus_and_ref_genome:
        in_paths = in_paths_by_gene_locus_and_ref_genome[(gene, locus, ref_genome)]
        out_paths = out_paths_by_gene_locus_and_ref_genome.get((gene, locus, ref_genome), [])
        for i in in_paths:
            for o in out_paths:
                in_out_match = get_in_out_match(i, o, genes_lengths[gene], minimal_gap_ratio, maximal_gap_ratio)
                if in_out_match is not None:
                    matches_per_gene_and_ref_genome[(gene, ref_genome)].append(in_out_match)
    matches_per_gene_and_ref_genome = {gene_ref_genome: keep_best_matches(matches, iou_th=iou_th)
                                       for gene_ref_genome, matches in matches_per_gene_and_ref_genome.items()}
    log.info(
        f"found {sum((len(m) for m in matches_per_gene_and_ref_genome.values()))} matches for {len(matches_per_gene_and_ref_genome)} gene-reference-genome pairs")
    return matches_per_gene_and_ref_genome

def keep_best_matches(matches: List, sorting_func=lambda x: (-x.score, -(x.end - x.start), x.start),
                      iou_th=IN_OUT_MATCH_IOU_TH) -> List:
    # TODO double check that I can get more than one match
    # sorted ascending on negated score and length, so the best comes first while ties still go to the
    # leftmost match - reverse=True would order those by descending start instead
    return sau.non_max_suppression_single_class(matches, sorting_func=sorting_func, iou_th=iou_th, reverse=False)


def get_ref_genome_species_dict_from_metadata_path(metadata_path):
    """Maps every reference genome in the metadata table to its species.

    Falls back to the species embedded in the Lineage column ('...;s__Escherichia coli') when there is
    no species column or a row leaves it empty, and to the first column when none is named Genome.

    Streamed rather than read with pandas: the reference table is tens of MB of columns this needs
    three of.
    """
    with open(metadata_path, 'r') as f:
        # QUOTE_NONE so that a '"' anywhere in a lineage stays part of the field instead of quoting it
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        if not reader.fieldnames:
            return {}
        reader.fieldnames = [name.strip().lower() for name in reader.fieldnames]
        genome_field = 'genome' if 'genome' in reader.fieldnames else reader.fieldnames[0]

        ref_genome_species_dict = {}
        for row in reader:
            genome = row[genome_field]
            if genome:
                ref_genome_species_dict[genome] = row.get('species') or species_from_lineage(row.get('lineage'))
    return ref_genome_species_dict


def read_contig_to_genome_map(contig_to_genome_path) -> dict:
    """contig id -> the reference genome it belongs to, as written when the reference DB was built.

    Missing or unreadable means "no map": callers fall back to deriving the genome from the contig
    name, which is what a reference database GInGeR did not build itself (--sample-specific-references)
    still relies on.
    """
    if not contig_to_genome_path or not os.path.exists(contig_to_genome_path):
        return {}
    with open(contig_to_genome_path) as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        if not reader.fieldnames:
            return {}
        reader.fieldnames = [name.strip().lower() for name in reader.fieldnames]
        if 'contig' not in reader.fieldnames or 'genome' not in reader.fieldnames:
            log.warning(f'{contig_to_genome_path} has no contig/genome columns - '
                        f'falling back to deriving genomes from contig names')
            return {}
        return {row['contig']: row['genome'] for row in reader if row.get('contig')}


def genome_from_contig_name(contig_name: str) -> str:
    """The genome a reference contig belongs to, guessed from its name.

    The fallback for a reference database GInGeR did not build and so has no contig->genome map for.
    Assumes the {genome_id}_{contig_num} naming that README documents for
    --sample-specific-references, and that UHGG's MGYG000260594_1 follows.
    """
    return re.split(r'[._]', contig_name)[0]


def build_contig_species_lookup(metadata_path, contig_to_genome_path=None) -> dict:
    """contig id -> (genome, species), the lookup a match on a reference contig is resolved with.

    Built by composing the contig->genome map written when the reference database was assembled with
    the genome->species mapping from the reference metadata. Without a map the genome is guessed from
    the contig name, and the lookup answers for any contig rather than a fixed set, so it is returned
    as a callable-backed dict-alike below.
    """
    genome_to_species = get_ref_genome_species_dict_from_metadata_path(metadata_path)
    contig_to_genome = read_contig_to_genome_map(contig_to_genome_path)
    if contig_to_genome:
        log.info(f'resolving reference contigs with a {len(contig_to_genome)} contig map')
    else:
        log.info('no reference contig->genome map - deriving genomes from reference contig names')
    return ContigSpeciesLookup(genome_to_species, contig_to_genome)


class ContigSpeciesLookup:
    """Resolves a reference contig to the (genome, species) it belongs to.

    Not a plain dict because without a contig->genome map the genome is derived from the contig name
    on demand, so the set of answerable contigs is not known up front.
    """

    def __init__(self, genome_to_species: dict, contig_to_genome: dict):
        self.genome_to_species = genome_to_species
        self.contig_to_genome = contig_to_genome

    def genome(self, contig_name: str) -> str:
        if self.contig_to_genome:
            return self.contig_to_genome.get(contig_name, '')
        return genome_from_contig_name(contig_name)

    def resolve(self, contig_name: str) -> tuple:
        """(genome, species) for a reference contig. species is 'unknown_{contig}' when the contig
        cannot be traced to a genome the metadata knows, which is the signal that the reference
        database and the metadata table do not describe the same genomes."""
        genome = self.genome(contig_name)
        species = self.genome_to_species.get(genome) if genome else None
        return genome, species or f'unknown_{contig_name}'


def species_from_lineage(lineage) -> str:
    """The species of a GTDB-style lineage - the 's__' rank of 'd__Bacteria;...;s__Escherichia coli'."""
    match = re.search(r'(?:^|;)s__([^;]+)', lineage) if lineage else None
    return match.group(1) if match else ''

@pu.step_timing
def process_in_and_out_paths_to_results(in_path_mapping_to_ref_genomes, out_path_mapping_to_ref_genomes, genes_lengths,
                                        paths_pident_filtering_th, minimal_gap_ratio,
                                        maximal_gap_ratio, metadata_path, contig_to_genome_path=None,
                                        contig_species_lookup=None):
    log.info('parsing the mapping of in and out paths')
    # the caller may have built the lookup already - it is also what the context level CSV is
    # written with, and reading the reference metadata twice is pure cost
    ref_species_dict = contig_species_lookup or build_contig_species_lookup(metadata_path, contig_to_genome_path)
    # keyed by (gene, gene copy, reference genome) - see read_and_filter_path_matches_per_gene
    in_paths_by_gene_locus_and_ref_genome = read_and_filter_path_matches_per_gene(
        mc.PathRefGenomeMatch, in_path_mapping_to_ref_genomes, paths_pident_filtering_th, ref_species_dict)
    out_paths_by_gene_locus_and_ref_genome = read_and_filter_path_matches_per_gene(
        mc.PathRefGenomeMatch, out_path_mapping_to_ref_genomes, paths_pident_filtering_th, ref_species_dict)
    if len(in_paths_by_gene_locus_and_ref_genome) == 0 or len(out_paths_by_gene_locus_and_ref_genome) == 0:
        log.info(f'GInGeR found {len(in_paths_by_gene_locus_and_ref_genome)} matches for incoming paths and '
                 f'{len(out_paths_by_gene_locus_and_ref_genome)} matches for outgoing paths. '
                 f'No results will be produced')
        return []
    log.info('generating in-out matches')
    return get_all_in_out_matches(in_paths_by_gene_locus_and_ref_genome, out_paths_by_gene_locus_and_ref_genome,
                                  genes_lengths, minimal_gap_ratio, maximal_gap_ratio)
