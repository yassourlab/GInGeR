import unittest
from shutil import rmtree
import pickle
import os
import ginger.verify_context_candidates as vcc

from tests import helper

TEST_FILES = helper.get_filedir()


class TestVerifyContextsCandidates(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.matches_no_dups_path = f'{TEST_FILES}/matches_no_overlaps.pkl'

    # def get_in_out_match(i, o, gene_length, minimal_gap_ratio, maximal_gap_ratio):
    #     for field in ['gene', 'nodes_list', 'bug', 'strand']:
    #         if getattr(i, field) != getattr(o, field):
    #             return None
    #     start, end = extract_start_and_end(i, o)
    #     start_end_diff = end - start
    #     gap_ratio = start_end_diff / gene_length
    #     score = (i.score * i.path_length + o.score * o.path_length) / (i.path_length + o.path_length)
    #     if minimal_gap_ratio < gap_ratio < maximal_gap_ratio:
    #         return mc.InOutPathsMatch(i, o, start, end, gap_ratio, score, gene_length)
    #     return None
    # def test_get_in_out_match(self):


    # def test_keep_best_matches(self):
    #     # with open(self.matches_no_dups_path, 'rb') as f:
    #     #     matches_before_dedupping = pickle.load(f)
    #     matches_after_dedupping = vcc.keep_best_matches(matches_before_dedupping)
    #     self.assertEqual(len(matches_after_dedupping), len(matches_before_dedupping),
    #                      'number of matches should stay the same')


if __name__ == '__main__':
    unittest.main()
