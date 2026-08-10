from pafpy import PafRecord
from collections import namedtuple
import re

# One copy of a gene in the assembly - where it sits on a contig, 0-based half-open. Identifies the
# thing a context was cut from, so it doubles as a grouping key and as a slice of the contig.
GeneLocus = namedtuple('GeneLocus', ['contig', 'start', 'end'])

# Where a graph node aligned to a gap-containing contig. node is the oriented short name of the node
# running in the contig's direction ('7285+'); origin is where that node starts in contig coordinates,
# which is before the alignment starts, since minimap2 clips an alignment's ends.
NodePlacement = namedtuple('NodePlacement', ['node', 'origin', 'score'])

CONTEXT_NAME_FIELDS = ['gene', 'contig', 'start', 'end', 'match_score', 'nodes', 'path', 'side']


class PathRefGenomeMatch:
    def __init__(self, paf_line: PafRecord, contigs_to_species: dict):
        # extract_contexts_candidates.context_name builds this, '|'-separated with the gene first.
        # Splitting from the right leaves the gene whatever '|' it contains - SARG's and CARD's gene
        # names have several - and no other field can contain one.
        fields = paf_line.qname.rsplit('|', len(CONTEXT_NAME_FIELDS) - 1)
        if len(fields) != len(CONTEXT_NAME_FIELDS):
            raise ValueError(f'context name {paf_line.qname!r} does not have the '
                             f'{"|".join(CONTEXT_NAME_FIELDS)} fields it is written with')
        gene, contig, start, end, match_score, nodes, path, side = fields

        self.query_name = paf_line.qname
        self.gene = gene
        # the copy of the gene this context was cut from. Contexts are only ever paired within one,
        # so that a pair describes a stretch of sequence that is really in the assembly
        self.locus = GeneLocus(contig, int(start), int(end))
        self.gene_match_score = float(match_score)
        self.nodes_list = nodes
        self.path = path
        self.side = side

        self.path_length = paf_line.qlen
        self.path_start = paf_line.qstart
        self.path_end = paf_line.qend

        self.strand = str(paf_line.strand)

        self.ref_genome = paf_line.tname
        # reference contigs are named {genome}_{contig}, and the species metadata is keyed by genome
        genome = re.split(r'[._]', self.ref_genome)[0]
        self.species = contigs_to_species.get(genome, f'unknown_{self.ref_genome}')
        self.ref_genome_length = paf_line.tlen
        self.ref_genome_start = paf_line.tstart
        self.ref_genome_end = paf_line.tend

        self.score = paf_line.mlen / self.path_length

    def __str__(self):
        return f'{self.ref_genome} {self.species} {self.gene} nodes: {self.nodes_list} path: {self.path} match score:{self.score} strand: {self.strand} {self.ref_genome_start} to {self.ref_genome_end}'


class GeneContigMatch:
    """A single alignment of a gene to a contig.

    mmseqs2 reports tstart/tend 1-based inclusive, with tstart > tend on the minus strand. Normalized
    here to the 0-based half-open convention every consumer slices with, so that
    contig_seq[match.start:match.end] is the gene and match.start is an offset into a node or segment.
    """

    def __init__(self, mmseq_line: str):
        target, query, tstart, tend, nident, qlen = mmseq_line.split('\t')
        self.gene_length = int(qlen) * 3
        self.gene = query.split(' ')[0]
        start_int = int(tstart)
        end_int = int(tend)

        self.strand = '+' if start_int < end_int else '-'

        self.contig = target
        self.start = min(start_int, end_int) - 1
        self.end = max(start_int, end_int)

        self.score = int(nident) / int(qlen)
        self.nodes_list = None
        self.start_in_first_node = None

    @property
    def aligned_length(self):
        """How much of the contig the gene actually covers - what anything looking for the end of the
        gene on the contig wants.

        NOT gene_length, which is the full reference protein (qlen * 3): mmseqs2 runs with -c 0.8, so an
        alignment may cover as little as 80% of it and may contain gaps. The two are equal only for a
        full-length ungapped hit.
        """
        return self.end - self.start

    def __str__(self):
        return f'{self.gene} {self.contig} {self.score} {self.nodes_list}'


class InOutPathsMatch:
    def __init__(self, in_path, out_path, start, end, gap_ratio, score, gene_length, gene=None, ref_genome=None, gene_match_score=None, in_context_score=None, out_context_score=None, locus=None):
        # the copy of the gene both contexts were cut from - they are only ever paired within one,
        # so the pair describes a real stretch of the assembly
        self.locus = locus
        self.in_path = in_path
        self.out_path = out_path
        if gene is None:
            self.gene = self.in_path.gene
        else:
            self.gene = gene

        if ref_genome is None:
            self.ref_genome = self.in_path.ref_genome
        else:
            self.ref_genome = ref_genome

        self.start = start
        self.end = end
        self.gap_ratio = gap_ratio
        self.score = score
        self.gene_length = gene_length
        self.gene_match_score = gene_match_score
        self.in_context_score = in_context_score
        self.out_context_score = out_context_score
