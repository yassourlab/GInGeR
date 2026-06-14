import unittest
import os
import pyfastg
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
        cls.context_level_output_path_with_plasmid_score = f'{TEST_FILES}/context_level_matches_with_plasmid_score.csv'
        cls.metadata_path = 'ginger/UHGG-metadata.tsv' # use this for running tests on github CI
        # cls.metadata_path = '../ginger/UHGG-metadata.tsv' # use this for running the test locally
        cls.species_level_output_path = 'test_species_level_matches.csv'
        cls.species_level_output_path_with_dups = f'test_species_level_matches_with_dups.csv'  # dups means that the same gene is matched to the same genome twice
        cls.species_level_output_path_with_plasmid_score = 'test_species_level_matches_with_plasmid_score.csv'
        # GT files
        cls.species_level_output_path_gt = f'{TEST_FILES}/species_level_matches.csv'
        cls.species_level_output_path_with_dups_gt = f'{TEST_FILES}/species_level_matches_with_dups.csv'

    @classmethod
    def tearDownClass(cls) -> None:
        if os.path.exists(cls.species_level_output_path):
            os.remove(cls.species_level_output_path)
        if os.path.exists(cls.species_level_output_path_with_plasmid_score):
            os.remove(cls.species_level_output_path_with_plasmid_score)

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

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv_with_plasmid_score(self):
        species_level_output = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_plasmid_score,
            self.metadata_path,
            self.species_level_output_path_with_plasmid_score, 1)

        self.assertTrue(os.path.exists(self.species_level_output_path_with_plasmid_score))
        self.assertIn('plasmid_score_mean', species_level_output.columns)
        self.assertIn('plasmid_score_most_common_context', species_level_output.columns)

        row = species_level_output.loc[('test_gene', 'Escherichia coli_D')]
        # plasmid_score_mean averages over the gene's two unique contexts (0.9 and 0.3)
        self.assertAlmostEqual(row['plasmid_score_mean'], 0.6)
        # plasmid_score_most_common_context is the score of the context matched to 2 genomes (0.9), not 1 (0.3)
        self.assertAlmostEqual(row['plasmid_score_most_common_context'], 0.9)

    def test_parse_paths_file(self):
        paths_file = f'{TEST_FILES}/SPAdes/contigs.paths'
        paths_w_gaps_file = f'{TEST_FILES}/SPAdes/contigs_w_added_gaps.paths'
        assembly_graph =  pyfastg.parse_fastg(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')

        parsed_paths_without_gaps, contigs_with_gaps = pu.parse_paths_file(paths_file, assembly_graph.nodes)
        self.assertEqual(len(parsed_paths_without_gaps), 2)
        self.assertEqual(len(contigs_with_gaps), 0)

        parsed_paths_without_gaps, contigs_with_gaps = pu.parse_paths_file(paths_w_gaps_file, assembly_graph.nodes)
        self.assertEqual(len(parsed_paths_without_gaps), 1)
        self.assertEqual(len(contigs_with_gaps), 1)


if __name__ == '__main__':
    unittest.main()
