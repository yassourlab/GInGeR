import os
import logging
import numpy as np
import pandas as pd
import timeit
from collections import defaultdict
from pafpy import PafFile

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


def parse_list_of_nodes(as_str):
    splt = as_str.split(',')
    return [node.replace(';', '') for node in splt]


def is_contig_name_func(line):
    return line.startswith('NODE')


def paf_record_to_dict(paf_record):
    # qstart matters: minimap2 clips the ends of an alignment, so the query does not necessarily
    # start where the alignment does, and anything placing the query in target coordinates has to
    # subtract it
    return dict(qname=paf_record.qname, qlen=paf_record.qlen, qstart=paf_record.qstart, qend=paf_record.qend,
                strand=str(paf_record.strand),
                tname=paf_record.tname, tlen=paf_record.tlen, tstart=paf_record.tstart, tend=paf_record.tend,
                mlen=paf_record.mlen)


def minimap_results_from_path(path, head_size=None):
    with open(path) as f:
        paf_file = PafFile(f)
        if head_size:
            head = [next(paf_file) for _ in range(head_size)]  # paf_file
        else:
            head = paf_file
        minimap_results = pd.DataFrame([paf_record_to_dict(paf_record) for paf_record in head])
    return minimap_results


def parse_path_segments(path_lines):
    """Splits the path lines of a single contig into segments, each a list of oriented graph nodes.
    SPAdes writes one line per segment, ending with ';' when another segment follows, but a segment
    may also be wrapped over several lines.
    """
    segments = []
    current_segment = []
    for line in path_lines:
        for part_index, part in enumerate(line.split(';')):
            if part_index and current_segment:  # the ';' preceding this part closed a segment
                segments.append(current_segment)
                current_segment = []
            if part.strip():
                current_segment += parse_list_of_nodes(part.strip())
    if current_segment:
        segments.append(current_segment)
    return segments


def parse_paths_file(paths_path, path_is_contig_name_func=is_contig_name_func):
    """Parses SPAdes' contigs.paths into {contig name: ordered list of path segments}.

    A contig assembled from a single graph path has a single segment. SPAdes splits a path with ';'
    when the contig was assembled from several graph paths joined using paired-end evidence - such a
    contig gets one segment per part and its name is also returned in the set of contigs with gaps
    (a join is real sequence in the contig, but it is not an edge of the graph, so the segments
    can't be stitched into one path).
    """
    contigs_to_path_lines = {}
    with open(paths_path) as paths_file:
        for line in paths_file.readlines():
            stripped_line = line.strip()
            if path_is_contig_name_func(stripped_line):
                contig_name = stripped_line
                contigs_to_path_lines[contig_name] = []
            else:
                contigs_to_path_lines[contig_name].append(stripped_line)

    parsed_paths = {contig_name: parse_path_segments(path_lines) for contig_name, path_lines in
                    contigs_to_path_lines.items()}
    contigs_with_gaps = {contig_name for contig_name, segments in parsed_paths.items() if len(segments) > 1}
    return parsed_paths, contigs_with_gaps


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


class PathGeometry:
    """Where each node of a graph path sits in the sequence that path spells out.

    Consecutive nodes of a path share a k-mer, so the path's sequence is its nodes concatenated with
    that overlap collapsed once per join. Everything that has to place something on a path is asking
    about the same offsets - where a gene starts, how long a contigs.paths segment is, how much of a
    context candidate the gene already covers - so they are computed once here.

    The overlaps are memoized because finding one is a search (get_sequence_overlap tries four
    likely k's and then scans), and the same path is walked again for every gene on the contig and
    for every context candidate enumerated off it.
    """

    def __init__(self, node_sequences):
        self._node_sequences = node_sequences
        self._overlaps = {}
        self._offsets = {}

    def _seq(self, node):
        return str(self._node_sequences[node].seq)

    def overlap(self, prev_node, node):
        if (prev_node, node) not in self._overlaps:
            try:
                self._overlaps[(prev_node, node)] = get_sequence_overlap(self._seq(prev_node), self._seq(node))
            except Exception as e:
                raise ValueError(f'no overlap between consecutive nodes {prev_node} and {node}: {e}') from e
        return self._overlaps[(prev_node, node)]

    def offsets(self, nodes):
        """Where each node starts in the path's sequence. The first is always 0."""
        if tuple(nodes) not in self._offsets:
            offsets = [0]
            for prev_node, node in zip(nodes, nodes[1:]):
                offsets.append(offsets[-1] + len(self._seq(prev_node)) - self.overlap(prev_node, node))
            self._offsets[tuple(nodes)] = offsets
        return self._offsets[tuple(nodes)]

    def length(self, nodes):
        return self.offsets(nodes)[-1] + len(self._seq(nodes[-1]))

    def sequence(self, nodes):
        seq = self._seq(nodes[0])
        for prev_node, node in zip(nodes, nodes[1:]):
            seq += self._seq(node)[self.overlap(prev_node, node):]
        return seq

    def nodes_covering(self, nodes, start, end):
        """The nodes of the path that [start, end) covers, and where start falls inside the first of
        them. Coordinates are the path's own.

        start_in_first_node is None when start lies outside the path, which leaves the gene
        unplaceable: an offset into the first node is what both context sides are measured from.
        """
        covering = []
        start_in_first_node = None
        for node, node_start in zip(nodes, self.offsets(nodes)):
            node_end = node_start + len(self._seq(node))
            if node_start < end and start < node_end:  # half-open, so a node the gene only abuts is not covered
                if start_in_first_node is None and node_start <= start:
                    start_in_first_node = start - node_start
                covering.append(node)
            elif covering:
                break
        return covering, start_in_first_node

    def holds_gene(self, nodes, start_in_first_node, aligned_length):
        """Whether a gene of aligned_length really sits at start_in_first_node on this path.

        The context extraction trims each side using these two numbers; when they don't hold it
        takes the wrong bases instead of failing, so they are checked before they are used.
        """
        if not nodes or start_in_first_node is None:
            return False
        gene_end = start_in_first_node + aligned_length
        return (0 <= start_in_first_node < len(self._seq(nodes[0]))  # the in side is measured inside the first node
                and gene_end <= self.length(nodes)  # the gene ends on the path, not past it
                and gene_end > self.offsets(nodes)[-1])  # and reaches the last node, which the out side is measured in


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





def compute_context_species_diversity(results_df, species_reference_counts,
                                       group_cols=('gene', 'in_context', 'out_context'),
                                       species_col='species', genome_col='Genome'):
    """For each unique in-gene-out trio (group_cols), compute the Shannon diversity index of the
    species it was matched to.

    For every species matched by a trio, the number of unique references it was matched to is
    divided by `species_reference_counts` (so species with more available references don't get
    more weight), and the corrected counts are normalized to probabilities before computing the
    Shannon diversity (-sum(p * ln(p))).
    """
    group_cols = list(group_cols)
    species_counts = results_df.groupby(group_cols + [species_col])[genome_col].nunique().rename(
        'n_genomes').reset_index()
    species_counts['corrected_count'] = species_counts['n_genomes'] / species_counts[species_col].map(
        species_reference_counts)
    total_corrected_count = species_counts.groupby(group_cols)['corrected_count'].sum().rename(
        'total_corrected_count')
    species_counts = species_counts.merge(total_corrected_count, on=group_cols)
    probabilities = species_counts['corrected_count'] / species_counts['total_corrected_count']
    species_counts['entropy_term'] = -probabilities * np.log(probabilities)
    # abs() avoids -0.0 (e.g. for trios matched to a single species, where the entropy term is -1*ln(1) = -0.0)
    diversity = species_counts.groupby(group_cols)['entropy_term'].sum().abs().rename(
        'context_species_diversity').reset_index()
    return results_df.merge(diversity, on=group_cols, how='left')


def _compute_single_context_confidence_score(results_df, species_reference_counts, context_col, species_col,
                                              score_col):
    """How confidently a single context column (in_context or out_context, considered on its
    own) points to each species it was matched to. See `compute_context_species_confidence_score`
    for the underlying logic.
    """
    species_counts = results_df.groupby([context_col, species_col]).size().rename('count').reset_index()
    species_counts['corrected_count'] = species_counts['count'] / species_counts[species_col].map(
        species_reference_counts)
    corrected_count_sum = species_counts.groupby(context_col)['corrected_count'].sum().rename('corrected_count_sum')
    species_counts = species_counts.merge(corrected_count_sum, on=context_col)
    species_counts[score_col] = species_counts['corrected_count'] / species_counts['corrected_count_sum']
    return species_counts[[context_col, species_col, score_col]]


def compute_context_species_confidence_score(results_df, species_reference_counts,
                                               context_cols=('in_context', 'out_context'),
                                               species_col='species'):
    """Compute how confidently a context (the in-gene-out trio) points to each species it was
    matched to, versus other species also matched by it.

    The in_context and out_context are scored independently using the same logic: for every
    species matched by the context, the number of matches is divided by
    `species_reference_counts` (so species with more available references don't get more
    weight), and the corrected counts are normalized so the scores of all species matched by
    that context sum to 1. `context_species_confidence_score` is the average of the in_context
    and out_context scores.
    """
    in_context_col, out_context_col = context_cols
    in_scores = _compute_single_context_confidence_score(results_df, species_reference_counts, in_context_col,
                                                          species_col, 'in_context_confidence_score')
    out_scores = _compute_single_context_confidence_score(results_df, species_reference_counts, out_context_col,
                                                           species_col, 'out_context_confidence_score')

    results_df = results_df.merge(in_scores, on=[in_context_col, species_col], how='left')
    results_df = results_df.merge(out_scores, on=[out_context_col, species_col], how='left')
    results_df['context_species_confidence_score'] = results_df[
        ['in_context_confidence_score', 'out_context_confidence_score']].mean(axis=1)
    return results_df.drop(columns=['in_context_confidence_score', 'out_context_confidence_score'])


@step_timing
def write_context_level_output_to_csv(output, csv_path: str, metadata_path: str, max_species_representatives: int):
    results_dict = defaultdict(list)
    for gene_species_tuple, matches_list in output.items():
        gene, reference = gene_species_tuple
        for match in matches_list:
            results_dict['gene'].append(gene)
            results_dict['reference_contig'].append(reference)
            results_dict['in_context'].append(match.in_path.query_name)
            results_dict['out_context'].append(match.out_path.query_name)
            # which copy of the gene in the assembly this pair of contexts flanks. both are cut from
            # the same one, so the row describes a stretch of sequence that is really on this contig
            results_dict['contig'].append(match.locus.contig if match.locus is not None else None)
            results_dict['gene_start_in_contig'].append(match.locus.start if match.locus is not None else None)
            results_dict['gene_end_in_contig'].append(match.locus.end if match.locus is not None else None)
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
    results_df = results_df.merge(metadata_df[metadata_cols_to_merge], on='Genome', how='left')

    if 'species' in results_df.columns:
        species_reference_counts = metadata_df['species'].value_counts().clip(upper=max_species_representatives)
        results_df = compute_context_species_diversity(results_df, species_reference_counts)
        results_df = compute_context_species_confidence_score(results_df, species_reference_counts)

    results_df.to_csv(csv_path, index=False)


@step_timing
def aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path, metadata_path,
                                                                         species_level_output_path,
                                                                         max_species_representatives: int,
                                                                         species_col='species'):
    context_level_df = pd.read_csv(context_level_output_path)
    metadata_df = pd.read_csv(metadata_path, sep='\t')

    context_level_df.columns = [x.lower() for x in context_level_df.columns]
    metadata_df.columns = [x.lower() for x in metadata_df.columns]

    group_cols = ['gene', species_col]
    context_cols = ['in_context', 'out_context']

    genomes_per_species = metadata_df[species_col].value_counts().to_frame()
    genomes_per_species[f'{species_col}_instances'] = genomes_per_species['count'].apply(
        lambda x: min(x, max_species_representatives))
    # combine in_context/out_context into a single column so their distinct-pair count can be
    # aggregated alongside genome/score in one groupby
    context_level_df['context_pair'] = list(zip(*(context_level_df[col] for col in context_cols)))
    agg_output = context_level_df.groupby(group_cols).aggregate(
        {'genome': ['nunique'], 'score': ['max'], 'context_pair': ['nunique']})
    agg_output.columns = ['_'.join(col) for col in agg_output.columns.values]
    agg_output = agg_output.merge(genomes_per_species, left_on=species_col, right_index=True, how='left')
    agg_output['references_ratio'] = agg_output['genome_nunique'] / agg_output[f'{species_col}_instances']
    agg_output = agg_output.rename(columns={'context_pair_nunique': 'n_contexts'})
    species_level_output = agg_output[['references_ratio', 'score_max', f'{species_col}_instances', 'n_contexts']]

    unique_contexts = context_level_df.drop_duplicates(group_cols + context_cols)

    if 'plasmid_score' in context_level_df.columns:
        # average plasmid score across the gene's unique contexts (don't over-weight contexts
        # that were matched to many reference genomes of the same species)
        plasmid_score_mean = unique_contexts.groupby(group_cols)['plasmid_score'].mean().rename('plasmid_score_mean')

        # plasmid score of the context(s) matched to the most reference genomes, averaging ties
        context_counts = context_level_df.groupby(group_cols + context_cols).size().reset_index(name='n_genomes')
        context_counts = context_counts.merge(unique_contexts[group_cols + context_cols + ['plasmid_score']],
                                               on=group_cols + context_cols)
        max_counts = context_counts.groupby(group_cols)['n_genomes'].transform('max')
        plasmid_score_most_common_context = context_counts[context_counts['n_genomes'] == max_counts].groupby(
            group_cols)['plasmid_score'].mean().rename('plasmid_score_most_common_context')

        species_level_output = species_level_output.join(plasmid_score_mean).join(plasmid_score_most_common_context)

    if 'context_species_confidence_score' in context_level_df.columns:
        # a context's confidence score is constant per (gene, species, in_context, out_context);
        # dedup before averaging so contexts matched to many reference genomes of the species
        # aren't over-weighted
        species_confidence_score = unique_contexts.groupby(group_cols)['context_species_confidence_score'].mean(
        ).rename('species_confidence_score')

        species_level_output = species_level_output.join(species_confidence_score)

    if species_level_output_path is not None:
        species_level_output.to_csv(species_level_output_path)
    return species_level_output
