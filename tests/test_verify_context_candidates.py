import os
import tempfile
import unittest

from pafpy import PafRecord

import ginger.verify_context_candidates as vcc
from ginger import matches_classes as mc

from tests import helper

TEST_FILES = helper.get_filedir()

FIRST_COPY = mc.GeneLocus('contig1', 100, 200)
SECOND_COPY = mc.GeneLocus('contig1', 8000, 8100)


def path_match(ref_genome_start, ref_genome_end, locus=FIRST_COPY, gene='geneA', ref_genome='ref1', strand='+',
               side='in'):
    """A context's alignment to a reference genome, as read out of a paf."""
    context_name = f'{gene}|{locus.contig}|{locus.start}|{locus.end}|1.0000|1+|1+|{side}'
    paf_line = (f'{context_name}\t100\t0\t100\t{strand}\t'
                f'{ref_genome}\t100000\t{ref_genome_start}\t{ref_genome_end}\t100\t100\t60')
    return mc.PathRefGenomeMatch(PafRecord.from_str(paf_line), {})


class GetAllInOutMatchesTest(unittest.TestCase):
    def test_contexts_are_only_paired_within_one_copy_of_the_gene(self):
        """Two copies of a gene, both matching the same reference genome close enough that pairing
        across them would pass the gap ratio filter.

        The incoming context of one copy followed by the outgoing context of the other describes a
        stretch of sequence that is on no contig, so it must not be reported.
        """
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100, FIRST_COPY)],
                    ('geneA', SECOND_COPY, 'ref1'): [path_match(9000, 9100, SECOND_COPY)]}
        out_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1200, 1300, FIRST_COPY, side='out')],
                     ('geneA', SECOND_COPY, 'ref1'): [path_match(9200, 9300, SECOND_COPY, side='out')]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        # the pairing is per copy, but the result is still keyed and deduplicated per gene and
        # reference genome, so both copies land under one key
        self.assertEqual(list(matches), [('geneA', 'ref1')])
        # without the cross pairing of the first copy's incoming context with the second's outgoing
        # one, which spans 1100-9200 and survives deduplication against both of these
        self.assertEqual([(match.start, match.end) for match in matches['geneA', 'ref1']],
                         [(1100, 1200), (9100, 9200)])

    def test_a_match_carries_the_copy_its_contexts_were_cut_from(self):
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100, FIRST_COPY)]}
        out_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1200, 1300, FIRST_COPY, side='out')]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        self.assertEqual([match.locus for match in matches['geneA', 'ref1']], [FIRST_COPY])

    def test_a_copy_with_no_outgoing_context_produces_no_match(self):
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100, FIRST_COPY)]}
        out_paths = {('geneA', SECOND_COPY, 'ref1'): [path_match(1200, 1300, SECOND_COPY, side='out')]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        self.assertEqual(matches, {})

    def test_contexts_of_different_copies_are_not_paired_even_within_one_group(self):
        # the grouping key is not the only thing keeping copies apart - the pair itself is checked,
        # so a mis-grouped context cannot slip through
        in_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1000, 1100, FIRST_COPY)]}
        out_paths = {('geneA', FIRST_COPY, 'ref1'): [path_match(1200, 1300, SECOND_COPY, side='out')]}

        matches = vcc.get_all_in_out_matches(in_paths, out_paths, {'geneA': 100}, 0, 100)

        self.assertEqual(matches, {})


class ReadAndFilterPathMatchesPerGeneTest(unittest.TestCase):
    IN_PATHS_PAF = f'{TEST_FILES}/in_paths_to_reference.paf'
    # every record in that paf is for this one context
    LOCUS = mc.GeneLocus('NODE_1_length_1000_cov_140.620106', 336, 615)

    def test_matches_are_grouped_by_the_copy_of_the_gene(self):
        grouped = vcc.read_and_filter_path_matches_per_gene(mc.PathRefGenomeMatch, self.IN_PATHS_PAF, 0.9, {})

        self.assertEqual(sorted(grouped), [('test_gene', self.LOCUS, 'MGYG000077121_281'),
                                           ('test_gene', self.LOCUS, 'MGYG000260594_1')])


class RefGenomeSpeciesDictTest(unittest.TestCase):
    """Every context level row's species comes from this dict, and the reference metadata does not
    always name the species in a column of its own."""

    def _dict_from(self, content):
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as f:
            f.write(content)
        try:
            return vcc.get_ref_genome_species_dict_from_metadata_path(f.name)
        finally:
            os.unlink(f.name)

    def test_reads_the_species_column(self):
        self.assertEqual(self._dict_from('Genome\tspecies\nG1\tE. coli\nG2\tS. aureus\n'),
                         {'G1': 'E. coli', 'G2': 'S. aureus'})

    def test_falls_back_to_the_lineage_when_the_species_is_empty(self):
        self.assertEqual(
            self._dict_from('Genome\tspecies\tLineage\nG1\t\td__Bacteria;g__Escherichia;s__Escherichia coli\n'),
            {'G1': 'Escherichia coli'})

    def test_falls_back_to_the_lineage_when_there_is_no_species_column(self):
        self.assertEqual(self._dict_from('Genome\tLineage\nG1\td__Bacteria;s__Blautia faecis\n'),
                         {'G1': 'Blautia faecis'})

    def test_a_lineage_with_no_species_rank_gives_no_species(self):
        self.assertEqual(self._dict_from('Genome\tspecies\tLineage\nG1\t\td__Bacteria;g__Escherichia\n'), {'G1': ''})

    def test_falls_back_to_the_first_column_when_none_is_named_genome(self):
        self.assertEqual(self._dict_from('Accession\tspecies\nA1\tE. coli\n'), {'A1': 'E. coli'})

    def test_headers_are_matched_ignoring_case_and_padding(self):
        self.assertEqual(self._dict_from(' GENOME \t SPECIES \nG1\tE. coli\n'), {'G1': 'E. coli'})

    def test_a_row_with_no_genome_is_skipped(self):
        self.assertEqual(self._dict_from('Genome\tspecies\nG1\tE. coli\n\tOrphan\n'), {'G1': 'E. coli'})

    def test_a_row_that_stops_short_of_the_species(self):
        self.assertEqual(self._dict_from('Genome\tspecies\tLineage\nG1\nG2\tE. coli\n'), {'G1': '', 'G2': 'E. coli'})

    def test_a_quote_in_a_field_is_not_treated_as_quoting(self):
        self.assertEqual(self._dict_from('Genome\tspecies\nG1\tstrain "x" sp.\n'), {'G1': 'strain "x" sp.'})

    def test_a_table_with_no_rows(self):
        self.assertEqual(self._dict_from('Genome\tspecies\n'), {})
        self.assertEqual(self._dict_from(''), {})


class KeepBestMatchesTest(unittest.TestCase):
    """keep_best_matches and the gene alignment NMS share one loop but break ties in opposite
    directions - so the direction is what has to be pinned."""

    class Match:
        def __init__(self, score, start, end, tag):
            self.score, self.start, self.end, self.tag = score, start, end, tag
            self.gene = tag

    def test_the_higher_scoring_of_two_overlapping_matches_wins(self):
        low = self.Match(0.90, 0, 100, 'low')
        high = self.Match(0.95, 5, 105, 'high')
        self.assertEqual([m.tag for m in vcc.keep_best_matches([low, high], iou_th=0.3)], ['high'])

    def test_matches_that_do_not_overlap_are_both_kept(self):
        left = self.Match(0.9, 0, 100, 'left')
        right = self.Match(0.9, 5000, 5100, 'right')
        kept = {m.tag for m in vcc.keep_best_matches([left, right], iou_th=0.3)}
        self.assertEqual(kept, {'left', 'right'})

    def test_an_exact_tie_goes_to_the_leftmost_match(self):
        # the gene alignment NMS sorts the same key with reverse=True and would keep 'right' instead
        left = self.Match(0.9, 0, 100, 'left')
        right = self.Match(0.9, 5, 105, 'right')
        for order in ([left, right], [right, left]):  # and the input's order must not decide it
            self.assertEqual([m.tag for m in vcc.keep_best_matches(order, iou_th=0.3)], ['left'])


if __name__ == '__main__':
    unittest.main()
