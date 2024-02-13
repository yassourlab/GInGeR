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

    # def test_keep_best_matches(self):
    #     # with open(self.matches_no_dups_path, 'rb') as f:
    #     #     matches_before_dedupping = pickle.load(f)
    #     matches_after_dedupping = vcc.keep_best_matches(matches_before_dedupping)
    #     self.assertEqual(len(matches_after_dedupping), len(matches_before_dedupping),
    #                      'number of matches should stay the same')


if __name__ == '__main__':
    unittest.main()
