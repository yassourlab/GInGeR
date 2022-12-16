import unittest
import os
from helper import get_test_data_dir,get_test_output_dir
from ginger import pipeline_utils as pu
from pathlib import Path

class PipelineUtilsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data_dir = get_test_data_dir()
        cls.context_level_output_path = f'{cls.data_dir}/context_level_matches.csv'
        cls.metadata_path = f'{Path(__file__).resolve().parent.parent}/ginger/UHGG-metadata.tsv'
        cls.species_level_output_path = f'{get_test_output_dir()}/test_species_level_matches.csv'

    @classmethod
    def tearDownClass(cls) -> None:
        os.remove(cls.species_level_output_path)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv(self):
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(self.context_level_output_path,
                                                                                self.metadata_path,
                                                                                self.species_level_output_path, 1)
        self.assertTrue(os.path.exists(self.species_level_output_path))

        with open(self.species_level_output_path, 'r') as f1, open(
                f'{self.data_dir}/species_level_matches.csv', 'r') as f2:
            lines_1 = f1.readlines()
            lines_2 = f2.readlines()
        self.assertListEqual(lines_1, lines_2)


if __name__ == '__main__':
    unittest.main()
