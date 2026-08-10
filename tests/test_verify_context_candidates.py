import unittest

from pafpy import PafRecord

import ginger.verify_context_candidates as vcc
from ginger import matches_classes as mc

from tests import helper

TEST_FILES = helper.get_filedir()

FIRST_COPY = mc.GeneLocus('contig1', 100, 200)
SECOND_COPY = mc.GeneLocus('contig1', 8000, 8100)


def path_match(ref_genome_start, ref_genome_end, gene='geneA', ref_genome='ref1', strand='+'):
    """A context's alignment to a reference genome, as read out of a paf."""
    paf_line = (f'{gene}_nodes_1+_match_1.0000_path_1+\t100\t0\t100\t{strand}\t'
                f'{ref_genome}\t100000\t{ref_genome_start}\t{ref_genome_end}\t100\t100\t60')
    return mc.PathRefGenomeMatch(PafRecord.from_str(paf_line), {})


class GetAllInOutMatchesTest(unittest.TestCase):
    def test_contexts_are_only_paired_within_one_copy_of_the_gene(self):
        """Two copies of a gene, both matching the same reference genome close enough that pairing
        across them would pass the gap ratio filter.

        The incoming context of one copy followed by the outgoing context of the other describes a
        stretch of sequence that is on no contig, so it must not be reported.
        """
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100)],
                    ('geneA', SECOND_COPY, 'ref1'): [path_match(9000, 9100)]}
        out_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1200, 1300)],
                     ('geneA', SECOND_COPY, 'ref1'): [path_match(9200, 9300)]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        # the pairing is per copy, but the result is still keyed and deduplicated per gene and
        # reference genome, so both copies land under one key
        self.assertEqual(list(matches), [('geneA', 'ref1')])
        # without the cross pairing of the first copy's incoming context with the second's outgoing
        # one, which spans 1100-9200 and survives deduplication against both of these
        self.assertEqual([(match.start, match.end) for match in matches['geneA', 'ref1']],
                         [(1100, 1200), (9100, 9200)])

    def test_a_match_carries_the_copy_its_contexts_were_cut_from(self):
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100)]}
        out_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1200, 1300)]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        self.assertEqual([match.locus for match in matches['geneA', 'ref1']], [FIRST_COPY])

    def test_a_copy_with_no_outgoing_context_produces_no_match(self):
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100)]}
        out_paths = {('geneA', SECOND_COPY, 'ref1'): [path_match(1200, 1300)]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        self.assertEqual(matches, {})


class ReadAndFilterPathMatchesPerGeneTest(unittest.TestCase):
    IN_PATHS_PAF = f'{TEST_FILES}/in_paths_to_reference.paf'
    # every record in that paf is for this one context
    CONTEXT_NAME = 'test_gene_nodes_5+_path_5+'
    LOCUS = mc.GeneLocus('NODE_1_length_1000_cov_140.620106', 336, 615)

    def test_matches_are_grouped_by_the_copy_of_the_gene(self):
        grouped = vcc.read_and_filter_path_matches_per_gene(mc.PathRefGenomeMatch, self.IN_PATHS_PAF, 0.9, {},
                                                            {self.CONTEXT_NAME: self.LOCUS})

        self.assertEqual(sorted(grouped), [('test_gene', self.LOCUS, 'MGYG000077121_281'),
                                           ('test_gene', self.LOCUS, 'MGYG000260594_1')])

    def test_matches_are_dropped_when_the_context_is_not_in_the_loci_table(self):
        # only possible if the paf and the table came from different runs. dropping the matches
        # beats guessing which copy of the gene they belong to
        grouped = vcc.read_and_filter_path_matches_per_gene(mc.PathRefGenomeMatch, self.IN_PATHS_PAF, 0.9, {}, {})

        self.assertEqual(grouped, {})


if __name__ == '__main__':
    unittest.main()
