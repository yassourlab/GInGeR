import unittest
import os

from ginger import pipeline_utils as pu

from tests import helper

TEST_FILES = helper.get_filedir()

"""
This test does not pass locally, but should pass on githubs CI.
When running it locally you will get the following error:
"FileNotFoundError: [Errno 2] No such file or directory: 'ginger/UHGG-metadata.tsv'"
If you do wish for it to pass locally (is a reasonable request), replace "cls.metadata_path = 'ginger/UHGG-metadata.tsv'"
with cls.metadata_path = '../ginger/UHGG-metadata.tsv' 
"""


class PipelineUtilsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.context_level_output_path = f'{TEST_FILES}/context_level_matches.csv'
        cls.context_level_output_path_with_dups = f'{TEST_FILES}/context_level_matches_with_dups.csv'  # dups means that the same gene is matched to the same genome twice
        cls.metadata_path = 'ginger/UHGG-metadata.tsv' # use this for running tests on github CI
        # cls.metadata_path = '../ginger/UHGG-metadata.tsv' # use this for running the test locally
        cls.species_level_output_path = 'test_species_level_matches.csv'
        cls.species_level_output_path_with_dups = f'test_species_level_matches_with_dups.csv'  # dups means that the same gene is matched to the same genome twice
        # GT files
        cls.species_level_output_path_gt = f'{TEST_FILES}/species_level_matches.csv'
        cls.species_level_output_path_with_dups_gt = f'{TEST_FILES}/species_level_matches_with_dups.csv'

    @classmethod
    def tearDownClass(cls) -> None:
        if os.path.exists(cls.species_level_output_path):
            os.remove(cls.species_level_output_path)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv(self):
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(self.context_level_output_path,
                                                                                self.metadata_path,
                                                                                self.species_level_output_path, 1)
        self.assertTrue(os.path.exists(self.species_level_output_path))

        with open(self.species_level_output_path, 'r') as out_f, open(self.species_level_output_path_gt, 'r') as gt_f:
            out_lines = out_f.readlines()
            gt_lines = gt_f.readlines()
        self.assertListEqual(out_lines, gt_lines)

        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_dups,
            self.metadata_path,
            self.species_level_output_path_with_dups, 1)

        with open(self.species_level_output_path_with_dups, 'r') as out_f, open(
                self.species_level_output_path_with_dups_gt, 'r') as gt_f:
            out_lines = out_f.readlines()
            gt_lines = gt_f.readlines()
        self.assertListEqual(out_lines, gt_lines)


if __name__ == '__main__':
    unittest.main()
